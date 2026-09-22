"""
Production-grade Site Characterization Planner & Geotechnical Boring Optimizer
Based on Canadian Foundation Engineering Manual (CFEM), NAS (2019), and Robertson (2006).
"""

import io
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import brentq
import streamlit as st

# -----------------------------------------------------------------------------
# 1. CORE THEORETICAL & COMPUTATIONAL ENGINE
# -----------------------------------------------------------------------------

@dataclass
class SoilLayer:
    name: str
    thickness: float
    gamma_bulk: float
    gamma_sat: float
    phi: float
    c: float
    is_sand: bool


class BoussinesqStressEngine:
    """Calculates 3D stress field under rectangular footing and overburden stresses."""

    @staticmethod
    def calculate_corner_stress_fraction(b_prime: float, l_prime: float, z: float) -> float:
        """
        Calculates stress influence factor under corner of rectangular area b_prime x l_prime.
        Exact Fadum / Boussinesq integration formulation.
        """
        if z <= 1e-4:
            return 0.25
        m = b_prime / z
        n = l_prime / z
        m2 = m * m
        n2 = n * n
        r2 = m2 + n2 + 1.0
        v1 = (m * n) / np.sqrt(r2)
        v2 = (m2 + n2 + 2.0) / ((m2 + 1.0) * (n2 + 1.0))
        term1 = v1 * v2
        arg2 = (m * n) / np.sqrt(m2 * n2 + r2)
        arg2 = np.clip(arg2, -1.0, 1.0)
        term2 = np.arcsin(arg2)
        influence = (1.0 / (2.0 * np.pi)) * (term1 + term2)
        return float(np.clip(influence, 0.0, 0.25))

    @classmethod
    def calculate_center_stress(cls, width: float, length: float, q_applied: float, z: float) -> float:
        """Stress increase under center of footing (superposition of 4 equal quarters)."""
        b_prime = width / 2.0
        l_prime = length / 2.0
        influence_quarter = cls.calculate_corner_stress_fraction(b_prime, l_prime, z)
        return float(4.0 * influence_quarter * q_applied)

    @staticmethod
    def calculate_geostatic_effective_stress(
        depth: float,
        gw_table: float,
        gamma_bulk: float = 18.5,
        gamma_sat: float = 20.0,
        gamma_w: float = 9.81
    ) -> float:
        """Calculates vertical effective stress with single or two-zone water table."""
        if depth <= 0.0:
            return 1e-3
        if depth <= gw_table:
            return depth * gamma_bulk
        sigma_v_dry = gw_table * gamma_bulk
        sigma_v_sub = (depth - gw_table) * (gamma_sat - gamma_w)
        return max(sigma_v_dry + sigma_v_sub, 1e-3)


class SitePlannerEngine:
    """Core geotechnical optimization and recommendation engine based on CFEM and NAS 2019."""

    def __init__(
        self,
        struct_type: str,
        length: float,
        width: float,
        depth_or_height: float,
        q_applied: float,
        gw_table: float,
        bedrock_depth: float,
        variability: str,
        risk_level: str,
    ):
        self.struct_type = struct_type
        self.length = max(length, 0.1)
        self.width = max(width, 0.1)
        self.depth_or_height = max(depth_or_height, 0.0)
        self.q_applied = max(q_applied, 1.0)
        self.gw_table = gw_table
        self.bedrock_depth = bedrock_depth
        self.variability = variability
        self.risk_level = risk_level
        self.footing_area = self.length * self.width

    def compute_elastic_termination_depth(self) -> Tuple[float, float, float]:
        """
        Determines the analytical depth where:
        1. Delta sigma < 10% q_applied
        2. Delta sigma < 5% sigma'_v0
        Returns (depth_10pct, depth_5pct, governing_depth).
        """
        def f_10pct(z: float) -> float:
            delta_sigma = BoussinesqStressEngine.calculate_center_stress(
                self.width, self.length, self.q_applied, z
            )
            return delta_sigma - (0.10 * self.q_applied)

        def f_5pct(z: float) -> float:
            delta_sigma = BoussinesqStressEngine.calculate_center_stress(
                self.width, self.length, self.q_applied, z
            )
            sigma_v0 = BoussinesqStressEngine.calculate_geostatic_effective_stress(
                z, self.gw_table
            )
            return delta_sigma - (0.05 * sigma_v0)

        # Solve root for 10% criteria
        try:
            z_low, z_high = 0.1, 150.0
            if f_10pct(z_low) * f_10pct(z_high) <= 0:
                z_10 = float(brentq(f_10pct, z_low, z_high))
            else:
                z_10 = 2.0 * self.width
        except Exception:
            z_10 = 2.0 * self.width

        # Solve root for 5% effective stress criteria
        try:
            z_low, z_high = 0.1, 150.0
            if f_5pct(z_low) * f_5pct(z_high) <= 0:
                z_5 = float(brentq(f_5pct, z_low, z_high))
            else:
                z_5 = 2.5 * self.width
        except Exception:
            z_5 = 2.5 * self.width

        governing = min(z_10, z_5)
        return z_10, z_5, governing

    def evaluate_investigation_scope(self) -> Dict:
        """Applies Table 5.1 (NAS 2019) and CFEM criteria to obtain complete specs."""
        z_10, z_5, z_elastic = self.compute_elastic_termination_depth()
        min_depth = 0.0
        n_borings = 0
        spacing_desc = ""
        depth_rule_desc = ""

        if self.struct_type == "Shallow Foundation":
            # NAS 2019 Table 5.1: 2B for L <= 2B
            geometric_depth = 2.0 * self.width if self.length <= 2.0 * self.width else 1.5 * self.width
            min_depth = max(geometric_depth, z_elastic)
            depth_rule_desc = (
                f"Governed by max of geometric rule ({geometric_depth:.1f} m) and "
                f"elastic stress threshold ({z_elastic:.1f} m where net stress < 10% applied or < 5% effective)."
            )

            # Borehole counts per Table 5.1 and Section 5.3.3
            if self.footing_area < 250.0:
                n_borings = 3
                spacing_desc = "Minimum 3 borings for small footprints (<250 m²)."
            elif self.footing_area <= 1000.0:
                n_borings = 4 if self.variability == "Uniform" else 5
                spacing_desc = "4 to 5 borings (corners + center) for area between 250 and 1000 m²."
            else:
                base_count = 5
                additional = int(math.ceil((self.footing_area - 1000.0) / 400.0))
                n_borings = base_count + additional
                if self.variability != "Uniform":
                    n_borings = int(math.ceil(n_borings * 1.25))
                spacing_desc = "Grid pattern spaced 20 to 30 m across the footprint."

        elif self.struct_type == "Deep Foundation (Piles/Shafts)":
            pile_tip_depth = self.depth_or_height
            # NAS 2019: In soil: Extend below tip min 6m or 2x group dimension
            soil_extent = max(6.0, 2.0 * self.width)
            min_depth = pile_tip_depth + soil_extent
            depth_rule_desc = (
                f"Pile tip estimated at {pile_tip_depth:.1f} m. Investigation extends {soil_extent:.1f} m "
                f"below tip (max of 6 m and 2x group width {2.0 * self.width:.1f} m)."
            )

            if self.width <= 30.0:
                n_borings = 3
                spacing_desc = "3 locations per foundation element (<30 m width)."
            else:
                n_borings = 5
                spacing_desc = "4 to 5 locations per foundation element (>30 m width)."

            if self.variability != "Uniform":
                n_borings += 2
                spacing_desc += " Additional borings added due to variable subsurface conditions."

        elif self.struct_type == "Retaining Wall / Soil Nailing / Excavation":
            wall_height = self.depth_or_height
            min_depth = 2.0 * wall_height
            depth_rule_desc = (
                f"Investigation depth equals 2.0 x wall height ({2.0 * wall_height:.1f} m) below excavation base, "
                f"extending through all soft/compressible strata."
            )

            wall_len = max(self.length, 10.0)
            if wall_len <= 30.0:
                n_borings = 2
                spacing_desc = "Minimum 2 borings for wall lengths under 30 m."
            else:
                spacing = 35.0 if self.variability != "Uniform" else 55.0
                n_borings = max(2, int(math.ceil(wall_len / spacing)) + 1)
                spacing_desc = f"Boreholes spaced {spacing:.0f} m alternating along wall face and retention zone."

        elif self.struct_type == "Embankment Foundation":
            emb_height = self.depth_or_height
            min_depth = 2.0 * emb_height
            depth_rule_desc = (
                f"Investigation depth equals 2.0 x embankment height ({2.0 * emb_height:.1f} m), "
                f"extended if deep soft strata are present."
            )

            emb_length = max(self.length, 30.0)
            spacing = 60.0 if self.variability != "Uniform" else 120.0
            longitudinal_count = max(2, int(math.ceil(emb_length / spacing)) + 1)
            # 3 locations across critical transverse cross section
            n_borings = longitudinal_count + 2
            spacing_desc = f"Spaced {spacing:.0f} m along centerline plus 3 borings at critical transverse sections."

        # Bedrock coring rule check
        bedrock_rule_applied = False
        rock_core_length = 0.0
        if self.bedrock_depth < min_depth:
            min_depth = self.bedrock_depth + 3.0
            bedrock_rule_applied = True
            rock_core_length = 3.0
            depth_rule_desc += (
                f" Bedrock encountered at {self.bedrock_depth:.1f} m. Overburden drilled to rock and penetrated "
                f"minimum 3.0 m into sound bedrock to verify non-boulder stratigraphy."
            )

        # In-situ testing suite based on Robertson (2006) and CFEM Ch 5
        recommended_tests = []
        if self.risk_level == "Low Risk":
            recommended_tests = [
                "Continuous visual classification and stratigraphic soil logging.",
                "Standard Penetration Test (SPT) at 1.5 m depth intervals.",
                "Index testing: Water content, Grain size sieve analysis, and Atterberg limits.",
                "Pocket penetrometer or Torvane on cohesive split-spoon samples."
            ]
        elif self.risk_level == "Moderate Risk":
            recommended_tests = [
                "Electric Piezocone Penetrometer (CPTu) soundings with continuous porewater pressure (u2).",
                "SPT with calibrated energy ratio (ER_r) records and split-spoon sample recovery.",
                "Field Vane Shear Test (FVT) if soft cohesive deposits are encountered.",
                "Flat Dilatometer Test (DMT) for constrained deformation modulus and lateral stress (K0).",
                "Thin-walled open tube (Shelby tube, ASTM D1587) sampling in cohesive layers."
            ]
        else:  # High Risk
            recommended_tests = [
                "Seismic Piezocone (SCPTu) for continuous shear wave velocity (Vs) and small-strain G0.",
                "Self-Boring Pressuremeter (SBPMT) for in-situ horizontal stress (sigma_h0) and G-gamma decay.",
                "Porewater Pressure Dissipation (PPD) tests to evaluate horizontal consolidation (ch) and permeability (k).",
                "High-recovery rock coring (minimum NQ/HQ triple-tube wireline) with televiewer/GSI logging.",
                "Stationary Piston Sampler (ASTM D6519) or Block Sampling (Sherbrooke) for undisturbed Class 1 samples.",
                "Compulsory SPT energy calibration (ASTM D4633) prior to mobilization."
            ]

        return {
            "min_depth": round(min_depth, 2),
            "z_elastic_10pct": round(z_10, 2),
            "z_elastic_5pct": round(z_5, 2),
            "n_borings": int(n_borings),
            "spacing_description": spacing_desc,
            "depth_rationale": depth_rule_desc,
            "bedrock_rule_applied": bedrock_rule_applied,
            "rock_core_length": rock_core_length,
            "recommended_in_situ_tests": recommended_tests,
        }

    def generate_borehole_coordinates(self, n_borings: int) -> pd.DataFrame:
        """Generates realistic plan layout coordinates for boring schedule."""
        coords = []
        l = self.length
        w = self.width

        if self.struct_type in ["Shallow Foundation", "Deep Foundation (Piles/Shafts)"]:
            if n_borings <= 3:
                coords = [
                    {"Boring_ID": "BH-01", "X_coord": 0.2 * l, "Y_coord": 0.2 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-02", "X_coord": 0.8 * l, "Y_coord": 0.5 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-03", "X_coord": 0.3 * l, "Y_coord": 0.8 * w, "Elevation": 100.0},
                ]
            elif n_borings <= 5:
                coords = [
                    {"Boring_ID": "BH-01", "X_coord": 0.15 * l, "Y_coord": 0.15 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-02", "X_coord": 0.85 * l, "Y_coord": 0.15 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-03", "X_coord": 0.50 * l, "Y_coord": 0.50 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-04", "X_coord": 0.15 * l, "Y_coord": 0.85 * w, "Elevation": 100.0},
                    {"Boring_ID": "BH-05", "X_coord": 0.85 * l, "Y_coord": 0.85 * w, "Elevation": 100.0},
                ]
            else:
                nx = int(math.ceil(np.sqrt(n_borings * (l / w))))
                ny = int(math.ceil(n_borings / nx))
                xs = np.linspace(0.1 * l, 0.9 * l, nx)
                ys = np.linspace(0.1 * w, 0.9 * w, ny)
                idx = 1
                for x in xs:
                    for y in ys:
                        if idx <= n_borings:
                            coords.append(
                                {"Boring_ID": f"BH-{idx:02d}", "X_coord": round(x, 2), "Y_coord": round(y, 2), "Elevation": 100.0}
                            )
                            idx += 1
        elif self.struct_type == "Retaining Wall / Soil Nailing / Excavation":
            xs = np.linspace(0.0, l, n_borings)
            for i, x in enumerate(xs):
                y = 0.0 if (i % 2 == 0) else -1.2 * self.depth_or_height
                tag = "Wall-Face" if y == 0.0 else "Retained-Zone"
                coords.append(
                    {"Boring_ID": f"BH-{i+1:02d}", "X_coord": round(x, 2), "Y_coord": round(y, 2), "Elevation": 100.0, "Zone": tag}
                )
        else:  # Embankment
            xs = np.linspace(0.0, l, max(2, n_borings - 2))
            for i, x in enumerate(xs):
                coords.append({"Boring_ID": f"BH-{i+1:02d}", "X_coord": round(x, 2), "Y_coord": 0.0, "Elevation": 100.0, "Zone": "Centerline"})
            # Transverse points
            coords.append({"Boring_ID": f"BH-{len(coords)+1:02d}", "X_coord": round(0.5 * l, 2), "Y_coord": -0.8 * w, "Elevation": 100.0, "Zone": "Toe-Left"})
            coords.append({"Boring_ID": f"BH-{len(coords)+1:02d}", "X_coord": round(0.5 * l, 2), "Y_coord": 0.8 * w, "Elevation": 100.0, "Zone": "Toe-Right"})

        df = pd.DataFrame(coords[:n_borings])
        if "Zone" not in df.columns:
            df["Zone"] = "Footprint"
        return df


# -----------------------------------------------------------------------------
# 2. DATA SYNTHESIS & BENCHMARK PROFILE GENERATOR
# -----------------------------------------------------------------------------

def generate_benchmark_subsurface_data(target_depth: float, gw_table: float) -> pd.DataFrame:
    """Creates realistic, high-resolution continuous CPTu and SPT synthetic borehole log."""
    depths = np.linspace(0.5, max(target_depth + 3.0, 10.0), 120)
    data = []

    for z in depths:
        # Layer 1: Silty Sand Fill / Upper Sand (0 - 4m)
        if z <= 4.0:
            layer = "Dense Sand Fill"
            qc = 8.0 + 1.2 * z + np.random.normal(0, 0.4)
            fs = qc * 0.009 + np.random.normal(0, 0.01)
            u2 = 9.81 * max(0.0, z - gw_table)
            n60 = int(np.clip(qc * 4.2 + np.random.normal(0, 1), 6, 45))
        # Layer 2: Medium Soft Silty Clay (4 - 12m)
        elif z <= 12.0:
            layer = "Soft to Firm Marine Silty Clay"
            qc = 1.2 + 0.15 * (z - 4.0) + np.random.normal(0, 0.1)
            fs = qc * 0.035 + np.random.normal(0, 0.005)
            # Dynamic excess pore pressure during penetration
            u2 = (9.81 * max(0.0, z - gw_table)) + (22.0 * (z - 4.0)) + np.random.normal(0, 5.0)
            n60 = int(np.clip(qc * 2.8 + np.random.normal(0, 1), 2, 12))
        # Layer 3: Glacial Till / Dense Cohesionless Sand (>12m)
        else:
            layer = "Glacial Silt-Sand Till"
            qc = 14.0 + 0.8 * (z - 12.0) + np.random.normal(0, 0.8)
            fs = qc * 0.018 + np.random.normal(0, 0.03)
            u2 = 9.81 * max(0.0, z - gw_table)
            n60 = int(np.clip(qc * 3.5 + np.random.normal(0, 2), 25, 60))

        qt = qc + (1.0 - 0.75) * (u2 / 1000.0)  # a_net = 0.75
        rf = (fs / qt) * 100.0 if qt > 0 else 0.0

        data.append({
            "Depth_m": round(z, 2),
            "Stratum": layer,
            "qc_MPa": round(max(qc, 0.2), 2),
            "qt_MPa": round(max(qt, 0.2), 2),
            "fs_MPa": round(max(fs, 0.005), 4),
            "u2_kPa": round(max(u2, 0.0), 1),
            "Rf_pct": round(max(rf, 0.1), 2),
            "SPT_N60": n60
        })

    return pd.DataFrame(data)


# -----------------------------------------------------------------------------
# 3. ADVANCED VISUALIZATIONS (PLOTLY & MATPLOTLIB)
# -----------------------------------------------------------------------------

class Visualizer:
    """Generates interactive Plotly figures and print-ready 300-DPI Matplotlib plots."""

    @staticmethod
    def create_interactive_stress_and_subsurface_plot(
        depth_array: np.ndarray,
        delta_sigma: np.ndarray,
        sigma_v0: np.ndarray,
        q_app: float,
        governing_depth: float,
        df_log: pd.DataFrame
    ) -> go.Figure:
        """Interactive 3-panel synchronized geotechnical dashboard."""
        fig = make_subplots(
            rows=1, cols=3,
            shared_yaxes=True,
            horizontal_spacing=0.06,
            subplot_titles=(
                "Stress Dissipation & Depth Criteria",
                "Corrected Cone Resistance (qt) Profile",
                "Pore Pressure (u2) & SPT N60 Profile"
            )
        )

        # Subplot 1: Boussinesq Stresses
        fig.add_trace(
            go.Scatter(x=delta_sigma, y=depth_array, mode='lines', name='Δσ (Boussinesq)',
                       line=dict(color='#1f77b4', width=2.5),
                       hovertemplate='Depth: %{y:.2f} m<br>Δσ: %{x:.1f} kPa'),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=0.10 * q_app * np.ones_like(depth_array), y=depth_array, mode='lines',
                       name='10% q_applied', line=dict(color='#d62728', dash='dash', width=1.5)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=0.05 * sigma_v0, y=depth_array, mode='lines',
                       name="5% σ'v0", line=dict(color='#2ca02c', dash='dot', width=1.5)),
            row=1, col=1
        )
        fig.add_hline(y=governing_depth, line=dict(color='#000000', width=2, dash='dashdot'),
                      annotation_text=f"Target Depth: {governing_depth:.1f} m", annotation_position="top right", row=1, col=1)

        # Subplot 2: CPTu qt
        fig.add_trace(
            go.Scatter(x=df_log["qt_MPa"], y=df_log["Depth_m"], mode='lines+markers', name='CPTu qt (MPa)',
                       marker=dict(size=3, color='#ff7f0e'), line=dict(color='#ff7f0e', width=2),
                       hovertemplate='Depth: %{y:.2f} m<br>qt: %{x:.2f} MPa'),
            row=1, col=2
        )

        # Subplot 3: Pore pressure u2 & SPT
        fig.add_trace(
            go.Scatter(x=df_log["u2_kPa"], y=df_log["Depth_m"], mode='lines', name='u2 (kPa)',
                       line=dict(color='#9467bd', width=2),
                       hovertemplate='Depth: %{y:.2f} m<br>u2: %{x:.1f} kPa'),
            row=1, col=3
        )
        fig.add_trace(
            go.Scatter(x=df_log["SPT_N60"] * 10.0, y=df_log["Depth_m"], mode='markers', name='SPT N60 (x10)',
                       marker=dict(symbol='square', size=5, color='#8c564b'),
                       hovertemplate='Depth: %{y:.2f} m<br>N60: %{marker.size} blows'),
            row=1, col=3
        )

        fig.update_yaxes(autorange='reversed', title_text="Depth Below Surface (m)", row=1, col=1)
        fig.update_xaxes(title_text="Stress (kPa)", row=1, col=1)
        fig.update_xaxes(title_text="qt (MPa)", row=1, col=2)
        fig.update_xaxes(title_text="u2 (kPa) / SPT N60*10", row=1, col=3)

        fig.update_layout(
            height=650,
            template="plotly_white",
            legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5),
            margin=dict(l=60, r=40, t=60, b=80)
        )
        return fig

    @staticmethod
    def create_interactive_plan_view(df_borings: pd.DataFrame, length: float, width: float) -> go.Figure:
        """Interactive 2D spatial layout showing planned boreholes and footprint."""
        fig = go.Figure()

        # Footprint boundary
        x_box = [0, length, length, 0, 0]
        y_box = [0, 0, width, width, 0]
        fig.add_trace(
            go.Scatter(x=x_box, y=y_box, mode='lines', name='Structure Perimeter',
                       fill="toself", fillcolor='rgba(200, 220, 240, 0.3)',
                       line=dict(color='#3366cc', width=2))
        )

        # Borehole locations
        fig.add_trace(
            go.Scatter(
                x=df_borings["X_coord"],
                y=df_borings["Y_coord"],
                mode='markers+text',
                name='Proposed Boreholes',
                text=df_borings["Boring_ID"],
                textposition="top center",
                marker=dict(size=12, color='#b30000', symbol='diamond-cross', line=dict(color='black', width=1.5)),
                hovertemplate='<b>%{text}</b><br>X: %{x:.1f} m<br>Y: %{y:.1f} m'
            )
        )

        fig.update_layout(
            title="Spatial Layout of Exploration Locations (NAS 2019 / CFEM)",
            xaxis_title="Easting / X Distance (m)",
            yaxis_title="Northing / Y Distance (m)",
            template="plotly_white",
            height=500,
            yaxis=dict(scaleanchor="x", scaleratio=1),
            margin=dict(l=50, r=50, t=50, b=50)
        )
        return fig

    @staticmethod
    def generate_static_publication_figure(
        depth_array: np.ndarray,
        delta_sigma: np.ndarray,
        sigma_v0: np.ndarray,
        q_app: float,
        governing_depth: float,
        df_log: pd.DataFrame
    ) -> io.BytesIO:
        """Generates a multi-panel, 300 DPI, publication-ready vector graphic."""
        plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
        fig, axes = plt.subplots(1, 3, figsize=(13, 7), sharey=True, dpi=300)

        # Panel 1: Stress Dissipation
        axes[0].plot(delta_sigma, depth_array, color='#1f77b4', lw=2.2, label=r'$\Delta\sigma_z$ (Center)')
        axes[0].axvline(0.10 * q_app, color='#d62728', ls='--', lw=1.5, label=r'$10\%\ q_{applied}$')
        axes[0].plot(0.05 * sigma_v0, depth_array, color='#2ca02c', ls=':', lw=1.8, label=r"$5\%\ \sigma'_{v0}$")
        axes[0].axhline(governing_depth, color='black', ls='-.', lw=1.5, label=f'Design Target ({governing_depth:.1f} m)')
        axes[0].set_xlabel('Stress (kPa)', fontsize=11, fontweight='bold')
        axes[0].set_ylabel('Depth Below Surface (m)', fontsize=11, fontweight='bold')
        axes[0].set_title('(a) Boussinesq Stress Infiltration', fontsize=12, fontweight='bold')
        axes[0].invert_yaxis()
        axes[0].legend(loc='lower right', frameon=True, fontsize=9)
        axes[0].grid(True, which='both', ls='--', alpha=0.6)

        # Panel 2: Corrected Cone Resistance
        axes[1].plot(df_log['qt_MPa'], df_log['Depth_m'], color='#ff7f0e', lw=2.0)
        axes[1].set_xlabel('Corrected Cone Tip Resistance $q_t$ (MPa)', fontsize=11, fontweight='bold')
        axes[1].set_title('(b) Synthetic Piezocone Profile', fontsize=12, fontweight='bold')
        axes[1].grid(True, which='both', ls='--', alpha=0.6)

        # Panel 3: Pore Pressure & SPT
        ax3_twin = axes[2].twiny()
        axes[2].plot(df_log['u2_kPa'], df_log['Depth_m'], color='#9467bd', lw=1.8, label='Pore Pressure $u_2$')
        ax3_twin.scatter(df_log['SPT_N60'], df_log['Depth_m'], color='#8c564b', s=25, marker='s', label='SPT $N_{60}$')
        axes[2].set_xlabel('Pore Pressure $u_2$ (kPa)', color='#9467bd', fontsize=11, fontweight='bold')
        ax3_twin.set_xlabel('SPT $N_{60}$ (Blows/0.3 m)', color='#8c564b', fontsize=11, fontweight='bold')
        axes[2].set_title('(c) In-Situ Water & Penetration Resistance', fontsize=12, fontweight='bold', pad=15)
        axes[2].grid(True, which='both', ls='--', alpha=0.6)

        fig.suptitle('Canadian Foundation Engineering Manual (CFEM) Site Characterization Suite', fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout()

        img_buf = io.BytesIO()
        plt.savefig(img_buf, format='png', dpi=300, bbox_inches='tight')
        img_buf.seek(0)
        plt.close(fig)
        return img_buf


# -----------------------------------------------------------------------------
# 4. EXCEL EXPORT ENGINE (OPENPYXL)
# -----------------------------------------------------------------------------

class ExcelReportExporter:
    """Generates an executive-level, professionally styled Excel (.xlsx) schedule."""

    @staticmethod
    def export_to_excel(
        plan_results: Dict,
        df_borings: pd.DataFrame,
        df_log: pd.DataFrame,
        input_params: Dict
    ) -> io.BytesIO:
        wb = openpyxl.Workbook()

        # Styles
        header_fill = PatternFill(start_color="1F497D", end_color="1F497D", fill_type="solid")
        section_fill = PatternFill(start_color="DCE6F1", end_color="DCE6F1", fill_type="solid")
        font_header = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
        font_section = Font(name="Calibri", size=11, bold=True, color="1F497D")
        font_title = Font(name="Calibri", size=14, bold=True, color="1F497D")
        font_data = Font(name="Calibri", size=10)
        thin_border = Border(
            left=Side(style='thin', color='B0B0B0'),
            right=Side(style='thin', color='B0B0B0'),
            top=Side(style='thin', color='B0B0B0'),
            bottom=Side(style='thin', color='B0B0B0')
        )

        # Sheet 1: Executive Summary & Scope
        ws_sum = wb.active
        ws_sum.title = "Investigation_Plan"
        ws_sum.views.sheetView[0].showGridLines = True

        ws_sum["A1"] = "SITE CHARACTERIZATION & DRILLING SPECIFICATION"
        ws_sum["A1"].font = font_title
        ws_sum["A2"] = "Compliant with Canadian Foundation Engineering Manual (CFEM) & NAS (2019)"
        ws_sum["A2"].font = Font(name="Calibri", size=10, italic=True)

        ws_sum["A4"] = "1. DESIGN INPUT PARAMETERS"
        ws_sum["A4"].font = font_section
        ws_sum["A4"].fill = section_fill

        inputs_data = [
            ("Structure / Element Type", input_params["struct_type"]),
            ("Footprint Dimensions (L x B)", f"{input_params['length']:.1f} m x {input_params['width']:.1f} m"),
            ("Foundation / Wall Depth or Height", f"{input_params['depth_or_height']:.1f} m"),
            ("Applied Net Foundation Stress (q0)", f"{input_params['q_applied']:.1f} kPa"),
            ("Groundwater Table Elevation / Depth", f"{input_params['gw_table']:.1f} m"),
            ("Estimated Depth to Sound Bedrock", f"{input_params['bedrock_depth']:.1f} m"),
            ("Subsurface Geological Variability", input_params["variability"]),
            ("Geotechnical Consequence Risk Class", input_params["risk_level"]),
        ]

        row = 5
        for param, val in inputs_data:
            ws_sum[f"A{row}"] = param
            ws_sum[f"B{row}"] = val
            ws_sum[f"A{row}"].font = font_data
            ws_sum[f"B{row}"].font = font_data
            ws_sum[f"A{row}"].border = thin_border
            ws_sum[f"B{row}"].border = thin_border
            row += 1

        row += 1
        ws_sum[f"A{row}"] = "2. DRILLING & CHARACTERIZATION SCOPE"
        ws_sum[f"A{row}"].font = font_section
        ws_sum[f"A{row}"].fill = section_fill
        row += 1

        results_data = [
            ("Recommended Minimum Borehole Depth", f"{plan_results['min_depth']:.2f} m"),
            ("Governing Depth Rationale", plan_results['depth_rationale']),
            ("Minimum Required Borehole Count", plan_results['n_borings']),
            ("Borehole Spacing & Layout Rules", plan_results['spacing_description']),
            ("Elastic 10% Applied Stress Depth", f"{plan_results['z_elastic_10pct']:.2f} m"),
            ("Elastic 5% Effective Stress Depth", f"{plan_results['z_elastic_5pct']:.2f} m"),
            ("Bedrock 3.0 m Core Infiltration Triggered", "YES" if plan_results['bedrock_rule_applied'] else "NO"),
        ]

        for desc, val in results_data:
            ws_sum[f"A{row}"] = desc
            ws_sum[f"B{row}"] = val
            ws_sum[f"A{row}"].font = font_data
            ws_sum[f"B{row}"].font = font_data
            ws_sum[f"A{row}"].border = thin_border
            ws_sum[f"B{row}"].border = thin_border
            row += 1

        row += 1
        ws_sum[f"A{row}"] = "3. MANDATED IN-SITU & FIELD TESTING SUITE"
        ws_sum[f"A{row}"].font = font_section
        ws_sum[f"A{row}"].fill = section_fill
        row += 1

        for i, test in enumerate(plan_results['recommended_in_situ_tests'], start=1):
            ws_sum[f"A{row}"] = f"Test Requirement {i}"
            ws_sum[f"B{row}"] = test
            ws_sum[f"A{row}"].font = font_data
            ws_sum[f"B{row}"].font = font_data
            ws_sum[f"A{row}"].border = thin_border
            ws_sum[f"B{row}"].border = thin_border
            row += 1

        # Sheet 2: Borehole Layout Schedule
        ws_bh = wb.create_sheet(title="Borehole_Schedule")
        ws_bh.views.sheetView[0].showGridLines = True
        bh_cols = ["Boring ID", "Zone / Element", "X Coordinate (m)", "Y Coordinate (m)", "Surface Elevation (m)", "Target Depth (m)"]
        ws_bh.append(bh_cols)
        for col_idx in range(1, len(bh_cols) + 1):
            c = ws_bh.cell(row=1, column=col_idx)
            c.fill = header_fill
            c.font = font_header
            c.alignment = Alignment(horizontal="center")

        for _, r_data in df_borings.iterrows():
            ws_bh.append([
                r_data["Boring_ID"],
                r_data.get("Zone", "Footprint"),
                r_data["X_coord"],
                r_data["Y_coord"],
                r_data["Elevation"],
                plan_results["min_depth"]
            ])

        for r in range(2, len(df_borings) + 2):
            for c in range(1, len(bh_cols) + 1):
                cell = ws_bh.cell(row=r, column=c)
                cell.font = font_data
                cell.border = thin_border
                cell.alignment = Alignment(horizontal="center")

        # Sheet 3: Subsurface Synthetic Log Data
        ws_log = wb.create_sheet(title="Synthetic_CPTu_Log")
        ws_log.views.sheetView[0].showGridLines = True
        log_cols = list(df_log.columns)
        ws_log.append(log_cols)
        for col_idx in range(1, len(log_cols) + 1):
            c = ws_log.cell(row=1, column=col_idx)
            c.fill = header_fill
            c.font = font_header
            c.alignment = Alignment(horizontal="center")

        for _, r_data in df_log.iterrows():
            ws_log.append(list(r_data))

        for r in range(2, len(df_log) + 2):
            for c in range(1, len(log_cols) + 1):
                cell = ws_log.cell(row=r, column=c)
                cell.font = font_data
                cell.border = thin_border

        # Adjust column widths
        for sheet in wb.worksheets:
            for col in sheet.columns:
                max_len = max(len(str(cell.value or '')) for cell in col)
                col_letter = get_column_letter(col[0].column)
                sheet.column_dimensions[col_letter].width = max(max_len + 3, 12)

        output_stream = io.BytesIO()
        wb.save(output_stream)
        output_stream.seek(0)
        return output_stream


# -----------------------------------------------------------------------------
# 5. STREAMLIT APPLICATION (ENTRY POINT)
# -----------------------------------------------------------------------------

def main():
    st.set_page_config(
        page_title="CFEM Geotechnical Site Characterization Planner",
        page_icon="🏗️",
        layout="wide"
    )

    st.title("🏗️ Geotechnical Site Characterization & Drilling Optimizer")
    st.markdown(
        """
        **Automated Exploration Scope Planner & In-Situ Testing Allocator**  
        *Strictly aligned with the Canadian Foundation Engineering Manual (CFEM), NAS (2019) Table 5.1, and Robertson (2006).*
        """
    )
    st.write("---")

    # Sidebar Engineering Inputs
    st.sidebar.header("📐 1. Structural & Geometric Inputs")
    struct_type = st.sidebar.selectbox(
        "Foundation / System Classification",
        [
            "Shallow Foundation",
            "Deep Foundation (Piles/Shafts)",
            "Retaining Wall / Soil Nailing / Excavation",
            "Embankment Foundation"
        ]
    )

    col_geom1, col_geom2 = st.sidebar.columns(2)
    with col_geom1:
        footing_width = st.number_input("Width B (m)", min_value=1.0, max_value=200.0, value=15.0, step=1.0)
    with col_geom2:
        footing_length = st.number_input("Length L (m)", min_value=1.0, max_value=500.0, value=25.0, step=1.0)

    dim_label = "Wall Height H (m)" if "Retaining" in struct_type or "Embankment" in struct_type else "Pile / Depth Df (m)"
    depth_or_height = st.sidebar.number_input(dim_label, min_value=1.0, max_value=80.0, value=12.0, step=1.0)
    q_applied = st.sidebar.number_input("Applied Net Bearing Stress q₀ (kPa)", min_value=10.0, max_value=2000.0, value=220.0, step=10.0)

    st.sidebar.header("🌍 2. Subsurface & Risk Framework")
    gw_table = st.sidebar.number_input("Groundwater Table Depth (m)", min_value=0.0, max_value=50.0, value=3.0, step=0.5)
    bedrock_depth = st.sidebar.number_input("Anticipated Bedrock Depth (m)", min_value=2.0, max_value=150.0, value=25.0, step=1.0)

    variability = st.sidebar.selectbox(
        "Geological Stratification Variability",
        ["Uniform / Predictable", "Uncertain / Highly Variable"]
    )

    risk_level = st.sidebar.selectbox(
        "Robertson (2006) Project Consequence Risk",
        ["Low Risk", "Moderate Risk", "High Risk"]
    )

    # Initialize Engine
    engine = SitePlannerEngine(
        struct_type=struct_type,
        length=footing_length,
        width=footing_width,
        depth_or_height=depth_or_height,
        q_applied=q_applied,
        gw_table=gw_table,
        bedrock_depth=bedrock_depth,
        variability=variability,
        risk_level=risk_level
    )

    plan_results = engine.evaluate_investigation_scope()
    df_borings = engine.generate_borehole_coordinates(plan_results["n_borings"])

    # Compute continuous stress curve
    z_eval = np.linspace(0.1, max(plan_results["min_depth"] * 1.3, 10.0), 100)
    delta_sig = np.array([
        BoussinesqStressEngine.calculate_center_stress(footing_width, footing_length, q_applied, z)
        for z in z_eval
    ])
    sigma_v0 = np.array([
        BoussinesqStressEngine.calculate_geostatic_effective_stress(z, gw_table)
        for z in z_eval
    ])

    # Upload or Benchmark Data
    st.sidebar.header("📂 3. Subsurface Stratigraphy Data")
    uploaded_file = st.sidebar.file_uploader("Upload Stratigraphic Log (CSV or Excel)", type=["csv", "xlsx"])
    if uploaded_file is not None:
        try:
            if uploaded_file.name.endswith(".csv"):
                df_log = pd.read_csv(uploaded_file)
            else:
                df_log = pd.read_excel(uploaded_file)
            st.sidebar.success("Custom borehole log loaded successfully!")
        except Exception as e:
            st.sidebar.error(f"Error parsing file: {e}. Fallback to benchmark data.")
            df_log = generate_benchmark_subsurface_data(plan_results["min_depth"], gw_table)
    else:
        df_log = generate_benchmark_subsurface_data(plan_results["min_depth"], gw_table)

    # Input Parameters dictionary for exporting
    input_params_dict = {
        "struct_type": struct_type,
        "length": footing_length,
        "width": footing_width,
        "depth_or_height": depth_or_height,
        "q_applied": q_applied,
        "gw_table": gw_table,
        "bedrock_depth": bedrock_depth,
        "variability": variability,
        "risk_level": risk_level
    }

    # Summary Metrics
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Min. Exploration Depth", f"{plan_results['min_depth']:.1f} m", delta="Governing Threshold")
    with m2:
        st.metric("Total Borehole Count", f"{plan_results['n_borings']} Holes", delta=f"{plan_results['n_borings'] * plan_results['min_depth']:.0f} m Total Drilling")
    with m3:
        st.metric("Elastic Infiltration (10% q₀)", f"{plan_results['z_elastic_10pct']:.1f} m")
    with m4:
        st.metric("Effective Stress (5% σ'v0)", f"{plan_results['z_elastic_5pct']:.1f} m")

    st.write("")

    # Tabbed Interface
    tab_overview, tab_interactive, tab_layout, tab_export = st.tabs([
        "📋 Specification & Guidelines",
        "📊 Stress & Subsurface Profiles",
        "🗺️ Spatial Layout & Boring Schedule",
        "📥 Professional Export Suite"
    ])

    with tab_overview:
        col_ov1, col_ov2 = st.columns([3, 2])
        with col_ov1:
            st.subheader("Governing Criteria & Design Rationale")
            st.info(f"**Depth Determination Rationale:** {plan_results['depth_rationale']}")
            st.success(f"**Spacing & Arrangement Standard:** {plan_results['spacing_description']}")

            if plan_results["bedrock_rule_applied"]:
                st.warning(
                    f"⚠️ **Bedrock Coring Protocol Triggered:** Competent bedrock surface detected at {bedrock_depth:.1f} m. "
                    f"Per CFEM Section 5.3.2, each borehole must core at least 3.0 m into sound rock to confirm bedrock "
                    f"integrity and eliminate misinterpretation of floating boulders in glacial till."
                )

        with col_ov2:
            st.subheader(f"Mandated In-Situ Tests ({risk_level})")
            for t_item in plan_results["recommended_in_situ_tests"]:
                st.markdown(f"- ✅ **{t_item}**")

    with tab_interactive:
        st.subheader("Boussinesq Stress Infiltration vs. Piezocone Logging")
        fig_stress = Visualizer.create_interactive_stress_and_subsurface_plot(
            depth_array=z_eval,
            delta_sigma=delta_sig,
            sigma_v0=sigma_v0,
            q_app=q_applied,
            governing_depth=plan_results["min_depth"],
            df_log=df_log
        )
        st.plotly_chart(fig_stress, use_container_width=True)

    with tab_layout:
        st.subheader("Borehole Distribution Plan")
        fig_plan = Visualizer.create_interactive_plan_view(df_borings, footing_length, footing_width)
        st.plotly_chart(fig_plan, use_container_width=True)

        st.markdown("#### Planned Exploration Location Coordinates")
        st.dataframe(df_borings, use_container_width=True)

    with tab_export:
        st.subheader("Executive Geotechnical Deliverables")
        st.markdown("Download print-ready vector figures and comprehensive Excel specifications for client reporting.")

        col_ex1, col_ex2 = st.columns(2)
        with col_ex1:
            st.markdown("##### 📄 Publication Figure (300 DPI Multi-Panel)")
            img_buffer = Visualizer.generate_static_publication_figure(
                depth_array=z_eval,
                delta_sigma=delta_sig,
                sigma_v0=sigma_v0,
                q_app=q_applied,
                governing_depth=plan_results["min_depth"],
                df_log=df_log
            )
            st.image(img_buffer, caption="Static High-Resolution Print Preview (CFEM Format)", use_container_width=True)
            st.download_button(
                label="⬇️ Download High-Res Report Figure (PNG)",
                data=img_buffer,
                file_name="CFEM_Site_Characterization_Profile.png",
                mime="image/png"
            )

        with col_ex2:
            st.markdown("##### 📊 Formatted Drilling Specification (Excel)")
            st.markdown(
                """
                Includes:
                - **Investigation_Plan:** Executive summary of geometry, depth rationale, and testing requirements.
                - **Borehole_Schedule:** Coordinates, elevations, and targeted drilling depths.
                - **Synthetic_CPTu_Log:** Numerical continuous profiling dataset.
                """
            )
            excel_data = ExcelReportExporter.export_to_excel(
                plan_results=plan_results,
                df_borings=df_borings,
                df_log=df_log,
                input_params=input_params_dict
            )
            st.download_button(
                label="⬇️ Download Complete Specification (.xlsx)",
                data=excel_data,
                file_name="CFEM_Site_Investigation_Program.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )


if __name__ == "__main__":
    main()
