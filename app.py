"""
Streamlit app: Radiator Tube Cooling Analysis
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import io
import json
from radiator import calculations as calc
from radiator import io as rio
from radiator import plots as rplots

# Default constants
DEFAULT_VOLUME_ML = 20.0
DEFAULT_DENSITY = 998.0
DEFAULT_CP = 4180.0
DEFAULT_LENGTH_MM = 100.0
DEFAULT_WIDTH_MM = 16.0
DEFAULT_HEIGHT_MM = 3.0
DEFAULT_AMBIENT_FALLBACK = 25.0

st.set_page_config(page_title="Radiator Tube Cooling Analysis", layout="wide")

# Title and subtitle
st.title("Radiator Tube Cooling Analysis")
st.markdown("*Heat Transfer | Heat Flux | Cooling Rates | Curve Fitting | Experiment Comparison*")

# Sidebar: experiment settings
st.sidebar.header("Experiment settings")

st.sidebar.subheader("Water properties")
volume_ml = st.sidebar.number_input("Total water volume (mL)", min_value=0.0, value=float(DEFAULT_VOLUME_ML), step=1.0)
density = st.sidebar.number_input("Water density (kg/m³)", value=float(DEFAULT_DENSITY))
mass_override = st.sidebar.number_input("Direct measured water mass (kg) — optional (overrides volume)", value=0.0, min_value=0.0, step=0.01)
cp = st.sidebar.number_input("Specific heat cp (J/(kg·K))", value=float(DEFAULT_CP))

st.sidebar.subheader("Ambient conditions")
use_measured_T0 = st.sidebar.checkbox("Use measured T0 where available", value=True)
ambient_fallback_const = st.sidebar.number_input("Fallback ambient T (°C) if T0 absent", value=float(DEFAULT_AMBIENT_FALLBACK))

st.sidebar.subheader("Tube geometry")
length_mm = st.sidebar.number_input("Tube length (mm)", value=float(DEFAULT_LENGTH_MM))
width_mm = st.sidebar.number_input("External width (mm)", value=float(DEFAULT_WIDTH_MM))
height_mm = st.sidebar.number_input("External height (mm)", value=float(DEFAULT_HEIGHT_MM))
area_method = st.sidebar.radio("Exposed area calculation", options=["Rectangular perimeter (calc)", "Manual area entry"], index=0)
if area_method == "Manual area entry":
    area_manual = st.sidebar.number_input("Exposed area (m²)", value=0.0, format="%.6f")
else:
    area_manual = None

st.sidebar.subheader("Air & fan")
air_velocity_meta = st.sidebar.number_input("Air velocity (m/s) — metadata", value=0.0)

st.sidebar.subheader("Analysis settings")
report_interval_s = st.sidebar.number_input("Report interval (s) — optional (0=use all readings)", value=0)
analysis_window_start = st.sidebar.number_input("Analysis window start elapsed_s (leave 0 to use start)", value=0)
analysis_window_end = st.sidebar.number_input("Analysis window end elapsed_s (0=use end)", value=0)

st.sidebar.write("---")
if st.sidebar.button("Create saved_runs directory"):
    os.makedirs("saved_runs", exist_ok=True)
    st.sidebar.success("saved_runs/ directory created (if not already present).")

# File upload
st.subheader("Upload CSV files")
uploaded_files = st.file_uploader("Upload one or multiple CSVs (first file used for single analysis)", type=['csv'], accept_multiple_files=True)

# Mode select
mode = st.selectbox("Mode", ["Single Experiment Analysis", "Experiment Comparison"])

if not uploaded_files:
    st.info("Upload at least one CSV file or use the included sample_experiment.csv from the project folder.")
    # show sample data instructions
    st.markdown("You can download the included sample CSV in the repository and upload it here.")
    st.stop()

# Helper: load each file and prepare
experiments = []
for f in uploaded_files:
    try:
        df_raw, mapping = rio.load_csv_file(f)
        df_prep, messages = rio.validate_and_prepare(df_raw, mapping, require_elapsed=True)
        experiments.append({'filename': getattr(f, "name", "uploaded"), 'raw': df_raw, 'mapping': mapping, 'df': df_prep, 'messages': messages})
    except Exception as e:
        st.error(f"Failed to read {getattr(f,'name','file')}: {e}")

# Show detected columns summary for first experiment
st.subheader("Detected columns and validation")
for exp in experiments:
    st.markdown(f"**{exp['filename']}**")
    for m in exp['messages']:
        if m.lower().startswith("missing"):
            st.warning(m)
        else:
            st.success(m)

# Single Experiment Analysis
if mode == "Single Experiment Analysis":
    exp = experiments[0]
    df = exp['df'].copy()
    st.markdown(f"### Analysis of: {exp['filename']}")
    # Ensure T2 and T3 present
    if 'T2_C' not in df.columns or 'T3_C' not in df.columns:
        st.error("T2 or T3 column missing; cannot perform cooling calculations.")
        st.stop()

    # Ambient selection: per-row T0 if available and chosen, otherwise use fallback
    if use_measured_T0 and 'T0_C' in df.columns:
        df['ambient_used_C'] = df['T0_C']
        ambient_used_descr = "Using measured per-row T0"
    else:
        df['ambient_used_C'] = ambient_fallback_const
        ambient_used_descr = f"Using fallback ambient temperature: {ambient_fallback_const} °C"
    st.write(ambient_used_descr)

    # Compute process chain
    df = calc.compute_tmean_and_gap(df)
    df = calc.compute_interval_deltas(df)
    # mass calculation
    if mass_override and mass_override > 0:
        mass_kg = mass_override
        mass_source = "direct mass override"
    else:
        mass_kg = calc.volume_ml_to_mass_kg(volume_ml, density)
        mass_source = f"calculated from {volume_ml} mL and density {density} kg/m³"
    df = calc.compute_interval_energy_and_cumulative(df, mass_kg=mass_kg, cp_j_per_kgk=cp)
    df = calc.compute_interval_heat_loss_rate(df)
    avg_net_W = calc.average_net_heat_loss_rate_over_window(df)
    # area calculation
    if area_manual and area_manual > 0:
        area_m2 = area_manual
        area_descr = f"Manual area provided: {area_m2:.6f} m²"
    else:
        area_m2 = calc.rectangular_tube_area_m2(length_mm, width_mm, height_mm)
        area_descr = f"Rectangular perimeter area: {area_m2:.6f} m² (2*(W+H)*L), mm->m conversions applied"
    # apparent flux
    apparent_flux = calc.compute_apparent_heat_flux(avg_net_W, area_m2) if avg_net_W is not None else None

    # Ambient-corrected
    df['ambient_used_C'] = df['ambient_used_C']
    df = calc.ambient_corrected_theta(df, ambient_col='ambient_used_C')

    # Fit models
    # Linear fit on Tmean vs time (s)
    linear_res = calc.linear_fit_time(df['Tmean'].values, df['elapsed_s'].values)
    # Exponential fit on theta
    exp_res = calc.exponential_fit_theta(df['theta_C'].values, df['elapsed_s'].values)

    # Top metric cards (display in columns)
    st.markdown("### Key metrics")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        total_energy_J = df['cumulative_energy_J'].iloc[-1] if 'cumulative_energy_J' in df.columns else np.nan
        st.metric("Sensible energy removed (J)", f"{float(total_energy_J):.1f}")
    with c2:
        st.metric("Average net heat-loss rate (W)", f"{float(avg_net_W):.3f}" if avg_net_W is not None else "—")
    with c3:
        st.metric("Exposed tube area (m²)", f"{area_m2:.6f}")
    with c4:
        st.metric("Apparent average heat flux (W/m²)", f"{float(apparent_flux):.3f}" if apparent_flux is not None else "—")

    # Calculation summary table
    st.markdown("### Calculation summary")
    summary = {
        "Analysis duration (s)": df['elapsed_s'].iloc[-1] - df['elapsed_s'].iloc[0],
        "Initial T2 (°C)": df['T2_C'].iloc[0],
        "Final T2 (°C)": df['T2_C'].iloc[-1],
        "Initial T3 (°C)": df['T3_C'].iloc[0],
        "Final T3 (°C)": df['T3_C'].iloc[-1],
        "Initial Tmean (°C)": df['Tmean'].iloc[0],
        "Final Tmean (°C)": df['Tmean'].iloc[-1],
        "Mean Tmean drop (°C)": df['Tmean'].iloc[0] - df['Tmean'].iloc[-1],
        "T2 avg cooling rate (°C/h)": df['rate_T2_C_per_h'].mean(skipna=True),
        "T3 avg cooling rate (°C/h)": df['rate_T3_C_per_h'].mean(skipna=True),
        "Tmean avg cooling rate (°C/h)": df['rate_Tmean_C_per_h'].mean(skipna=True),
        "Mean abs T2-T3 gap (°C)": df['Tgap_abs'].mean(skipna=True),
        "Max abs T2-T3 gap (°C)": df['Tgap_abs'].max(skipna=True),
        "Water mass used (kg)": mass_kg,
        "cp used (J/(kg·K))": cp,
        "Energy removed (J)": float(total_energy_J),
        "Average net heat-loss (W)": float(avg_net_W) if avg_net_W is not None else None,
        "Exposed area (m²)": area_m2,
        "Apparent average heat flux (W/m²)": apparent_flux
    }
    df_summary = pd.DataFrame.from_dict(summary, orient='index', columns=['value'])
    st.dataframe(df_summary)

    # Plots
    st.markdown("### Plots")
    colp1, colp2 = st.columns([2,1])
    with colp1:
        fig_temp = rplots.plot_temperature(df)
        st.plotly_chart(fig_temp, use_container_width=True)
        fig_theta = rplots.plot_theta(df)
        st.plotly_chart(fig_theta, use_container_width=True)
    with colp2:
        fig_energy = rplots.plot_energy_cumulative(df)
        st.plotly_chart(fig_energy, use_container_width=True)
        fig_rates = rplots.plot_cooling_rates(df)
        st.plotly_chart(fig_rates, use_container_width=True)

    # Fitted models display
    st.markdown("### Curve fits")
    if linear_res.get('success', True):
        st.subheader("Linear fit")
        st.write(f"Tmean = {linear_res['intercept']:.4f} + {linear_res['slope_per_min']:.6f} °C/min × t")
        st.write(f"Slope: {linear_res['slope_per_h']:.3f} °C/h; R² = {linear_res['R2']:.4f}; RMSE = {linear_res['RMSE_C']:.4f} °C")
    else:
        st.write("Linear fit not available:", linear_res.get('message'))

    if exp_res.get('success', False):
        st.subheader("Exponential fit (θ = A exp(-k t))")
        st.write(f"A = {exp_res['A']:.4f} °C; k = {exp_res['k_min_inv']:.6f} min⁻¹; τ = {exp_res['tau_min']:.2f} min")
        st.write(f"R² = {exp_res['R2']:.4f}; RMSE = {exp_res['RMSE_C']:.4f} °C")
        # overlay exponential onto theta plot
        xfit_min, yfit = exp_res['y_fit']
        st.plotly_chart(rplots.plot_theta(df).add_scatter(x=np.array(xfit_min)/1.0, y=yfit, mode='lines', name='Exponential fit'), use_container_width=True)
    else:
        st.write("Exponential fit not available or unsuccessful:", exp_res.get('message'))

    # Processed data table
    st.markdown("### Processed readings (expandable)")
    st.dataframe(df)

    # Download buttons
    st.markdown("### Downloads")
    csv_buf = df.to_csv(index=False).encode('utf-8')
    st.download_button("Download processed readings CSV", csv_buf, file_name=f"processed_{exp['filename']}", mime='text/csv')
    summary_csv = df_summary.reset_index().to_csv(index=False).encode('utf-8')
    st.download_button("Download summary CSV", summary_csv, file_name=f"summary_{exp['filename']}", mime='text/csv')

    # Save run to saved_runs/
    st.markdown("### Save analysis locally (saved_runs/)")
    run_label = st.text_input("Run label", value=f"run_{exp['filename'].replace('.','_')}")
    if st.button("Save run"):
        os.makedirs("saved_runs", exist_ok=True)
        # sanitize filename basics
        safe_label = "".join(c for c in run_label if c.isalnum() or c in (' ', '_', '-')).rstrip()
        base = f"saved_runs/{safe_label}"
        i = 0
        while os.path.exists(f"{base}_{i}_processed.csv"):
            i += 1
        processed_fn = f"{base}_{i}_processed.csv"
        summary_fn = f"{base}_{i}_summary.csv"
        meta_fn = f"{base}_{i}_meta.json"
        df.to_csv(processed_fn, index=False)
        df_summary.reset_index().to_csv(summary_fn, index=False)
        meta = {
            'label': run_label,
            'files_used': [exp['filename']],
            'settings': {
                'volume_ml': volume_ml, 'density': density, 'mass_override': mass_override, 'cp': cp,
                'length_mm': length_mm, 'width_mm': width_mm, 'height_mm': height_mm, 'area_m2': area_m2,
                'ambient_used_description': ambient_used_descr
            }
        }
        with open(meta_fn, 'w') as fh:
            json.dump(meta, fh, indent=2)
        st.success(f"Saved processed CSV, summary CSV and metadata JSON as {processed_fn}, {summary_fn}, {meta_fn}")
        st.warning("Saved files are in the app's local filesystem; storage on Streamlit Community Cloud may be ephemeral. Do not upload confidential data without checking policies.")
