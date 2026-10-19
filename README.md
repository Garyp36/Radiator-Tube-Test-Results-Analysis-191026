# Radiator Tube Cooling Analysis — Streamlit App

Overview
--------
This Streamlit application analyses radiator tube cooling experiments. It computes sensible energy removed from water, interval and average heat-loss rates, apparent heat flux (based on water sensible-energy change), cooling rates, ambient-corrected temperatures, and fits both linear and exponential cooling models.

IMPORTANT: T2 and T3 are endpoint temperatures. Their arithmetic mean (`Tmean = (T2 + T3)/2`) is an approximation and is NOT a true measured bulk temperature.

This project is intended for local use and research; do not deploy confidential data publicly without checking data policies (Streamlit Community Cloud storage may be ephemeral).

Files
-----
- app.py — Streamlit application entry point.
- radiator/calculations.py — calculation engine (verified with tests).
- radiator/io.py — CSV load and validation helpers.
- radiator/plots.py — Plotly plotting helpers.
- requirements.txt — Python dependencies.
- sample_experiment.csv — provided sample 90-minute experiment (10 readings).
- tests/test_calculations.py — automated tests for core calculation functions.
- .gitignore

Requirements
------------
- Python 3.9+ recommended (3.10 or 3.11 also fine).
- Pip to install dependencies.

Install
-------
1. Create virtual environment (recommended).

macOS / Linux:
```
python3 -m venv venv
source venv/bin/activate
```

Windows (PowerShell):
```
python -m venv venv
.\venv\Scripts\Activate.ps1
```

2. Install dependencies:
```
pip install -r requirements.txt
```

Run the app locally
-------------------
From the project root:
```
streamlit run app.py
```
Open the URL printed by Streamlit (typically http://localhost:8501).

CSV input format and accepted aliases
------------------------------------
Preferred columns:
`elapsed_s`, `datetime`, `T0_C`, `T1_C`, `T2_C`, `T3_C`, `target_C`, `heater_enabled`, `ssr_cmd`, `fan_pwm_pct`, `fan_rpm`, `air_speed_m_s`, `sensor_status`, `rtc_status`, `heater_fault`.

Accepted aliases:
- elapsed: `elapsed`, `seconds`, `time_s`, `timestamp` (or `datetime` to compute elapsed)
- T0: `T0`, `ambient`, `ambient_C`
- T1: `T1`, `heated_water`
- T2: `T2`, `inlet`, `inlet_C`
- T3: `T3`, `outlet`, `outlet_C`

Validation
----------
- The app detects which columns are present and reports missing required fields.
- `elapsed_s` (or `datetime`) and T2 & T3 are essential for temperature-based cooling analysis.
- Missing optional inputs disable the dependent calculations and the UI shows warnings.

Key calculations (implemented and used)
---------------------------------------
Units are indicated for each formula:

- Tmean (°C): (T2 + T3)/2 (labelled explicitly as mean of endpoints)
- Temperature gap (°C): |T2 - T3| (mean and maximum over window reported)
- Water mass (kg): if measured mass provided, it overrides; otherwise:
  mass = volume_mL * 1e-6 * density_kg_m3
- Interval cooling rate (°C/h): -ΔT / Δt_s × 3600
  (Negative ΔT for cooling -> positive cooling rate when temperature decreases)
- Interval energy removed (J): m * cp * (Tmean_prev - Tmean_current)
- Cumulative energy removed (J): m * cp * (Tmean_initial - Tmean_current)
- Interval average net heat-loss rate (W): interval_energy_J / dt_s
- Average net heat-loss rate over window (W): total_cumulative_energy / total_elapsed_duration
- Tube surface area (m²): rectangular exposed area = 2*(width + height)*length (all mm→m conversion)
- Apparent heat flux (W/m²): average_net_heat_loss_rate / exposed_area
  -- Labelled explicitly as "apparent/estimated" and energy-balance based; it is NOT a direct wall heat flux measurement.
- Ambient-corrected θ (°C): Tmean - T0 (per-row). Normalised θ = θ / θ_initial (guard against division by zero).

Curve fitting
-------------
- Linear fit: Tmean = intercept + slope × time (time in seconds). Reports slope in °C/min and °C/h, R² and RMSE.
- Exponential fit: θ(t) = A × exp(-k × t) with t in minutes. Uses scipy.optimize.curve_fit. Reports A, k (min⁻¹), τ=1/k (min), R², RMSE.

Limitations and assumptions
---------------------------
- Tmean is an approximation from two end sensors; it is not a direct bulk temperature measurement.
- Apparent heat flux is computed from the water's sensible-energy change only; it excludes heat storage in the aluminium tube, radiation and other heat paths.
- No flow sensor; water is assumed stationary.
- Do not extrapolate when computing threshold crossing times; interpolation between measured points is allowed.
- Normalised θ undefined if initial θ is zero or invalid — app returns blank and a warning in that case.

Experiment comparison
---------------------
- Upload at least two CSV files and select "Experiment Comparison".
- The app compares metrics and overlays plots. You can compare by:
  - Time to a common Tmean threshold (first downward crossing; linear interpolation between measured points).
  - Time to a common θ threshold.
  - Values at a common elapsed time (interpolation only if within the experiment range).
- If an experiment starts already below threshold, reported crossing time is 0. If never reached, report "Not reached".

Local storage
-------------
- Saving uses a `saved_runs/` directory created locally by the app.
- Saved outputs: processed CSV, summary CSV, meta JSON with settings and file names.
- Files are saved with unique names to avoid accidental overwrite.
- Streamlit Community Cloud storage may be ephemeral — the app warns about this.

Testing
-------
Run tests locally:
```
pytest -q
```
The tests included verify:
- volume-to-mass conversion
- rectangular area calculation
- energy and heat-loss computations
- cooling-rate sign
- linear and exponential fits on synthetic data
- threshold crossing

Packaging as a ZIP
------------------
From project root:
```
zip -r radiator_streamlit_app.zip . -x "venv/*" "*.pyc" "__pycache__/*" ".git/*"
```

Deployment to Streamlit Community Cloud
--------------------------------------
1. Create a GitHub repo with this project.
2. Push files to the repo.
3. On Streamlit Cloud, click "New app" → connect your GitHub repo → select branch and `app.py`.
4. Deploy. Note: saved_runs/ storage is ephemeral on Streamlit Cloud.

Privacy
-------
No authentication or cloud DB is implemented — do not use with confidential data in public deployments.
