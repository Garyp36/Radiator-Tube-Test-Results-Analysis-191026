"""
calculation functions for radiator_streamlit_app

Implements:
- column aliases mapping helpers
- mass calculation
- tmean and gap
- cooling rates
- energy removed (interval and cumulative)
- heat-loss rates
- area calculation
- apparent heat flux
- ambient-corrected temperatures
- linear and exponential fits with stats
"""

from typing import Optional, Tuple, Dict
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# -----------------------
# Column mapping helpers
# -----------------------
ALIASES = {
    'elapsed_s': ['elapsed_s', 'elapsed', 'seconds', 'time_s', 'sec'],
    'datetime': ['datetime', 'timestamp', 'time'],
    'T0_C': ['T0_C', 'T0', 'ambient_C', 'ambient', 'ambient_temp'],
    'T1_C': ['T1_C', 'T1', 'heated_water_C', 'heated_water'],
    'T2_C': ['T2_C', 'T2', 'inlet_C', 'inlet'],
    'T3_C': ['T3_C', 'T3', 'outlet_C', 'outlet'],
    'air_speed_m_s': ['air_speed_m_s', 'air_speed', 'air_m_s', 'air_mps'],
    'fan_pwm_pct': ['fan_pwm_pct', 'fan_pwm', 'fan_pwm_percent'],
    'fan_rpm': ['fan_rpm', 'rpm'],
}

def map_columns(df: pd.DataFrame) -> Dict[str,str]:
    """
    Detects columns from aliases and returns mapping standard_name -> actual_column
    """
    mapping = {}
    cols_lower = {c.lower(): c for c in df.columns}
    for standard, aliases in ALIASES.items():
        found = None
        for alias in aliases:
            if alias.lower() in cols_lower:
                found = cols_lower[alias.lower()]
                break
        if found:
            mapping[standard] = found
    # elapsed may be given by datetime only; caller handles that
    return mapping

# -----------------------
# Geometry & mass
# -----------------------

def volume_ml_to_mass_kg(volume_ml: float, density_kg_m3: float=998.0) -> float:
    """
    Convert volume in mL to mass in kg using density kg/m^3.

    mass (kg) = volume_ml * 1e-6 m^3 * density (kg/m^3)
    """
    return float(volume_ml) * 1e-6 * float(density_kg_m3)

def rectangular_tube_area_m2(length_mm: float, width_mm: float, height_mm: float) -> float:
    """Area (m²) = 2 * (width + height) * length; convert mm->m"""
    L = float(length_mm) / 1000.0
    W = float(width_mm) / 1000.0
    H = float(height_mm) / 1000.0
    return 2.0 * (W + H) * L

# -----------------------
# Data processing
# -----------------------

def prepare_dataframe(df: pd.DataFrame, mapping: Dict[str,str], require_elapsed: bool=True) -> pd.DataFrame:
    """
    Map aliases and create a working dataframe with at least:
    elapsed_s, T0_C (optional), T1_C (optional), T2_C, T3_C

    Sorting by elapsed_s, handling duplicate timestamps by keeping first occurrence (documented policy).
    """
    d = df.copy()
    # Column renames
    for std, col in mapping.items():
        if col != std:
            d = d.rename(columns={col: std})
    # If datetime exists but elapsed_s not present, try to compute elapsed_s from first timestamp
    if 'elapsed_s' not in d.columns and 'datetime' in d.columns:
        try:
            dt = pd.to_datetime(d['datetime'])
            base = dt.iloc[0]
            d['elapsed_s'] = (dt - base).dt.total_seconds()
        except Exception:
            # leave as-is; later validation will complain
            pass

    # Attempt to cast numeric columns
    for c in ['elapsed_s','T0_C','T1_C','T2_C','T3_C','air_speed_m_s']:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')

    # Drop rows that have neither elapsed nor datetime or missing both T2 and T3
    # We keep rows where at least one of T2 or T3 exists, but calculations need both.
    if require_elapsed:
        d = d[~d['elapsed_s'].isna()]
    # Sort by elapsed_s
    d = d.sort_values('elapsed_s').reset_index(drop=True)

    # Remove exact-duplicate elapsed_s keeping first occurrence (policy)
    if 'elapsed_s' in d.columns:
        d = d[~d['elapsed_s'].duplicated(keep='first')].reset_index(drop=True)

    return d

# -----------------------
# Core calculations
# -----------------------

def compute_tmean_and_gap(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds columns:
    - Tmean (mean of endpoints)
    - Tgap_abs = abs(T2 - T3)
    """
    d = df.copy()
    d['Tmean'] = np.nan
    if 'T2_C' in d.columns and 'T3_C' in d.columns:
        d['Tmean'] = (d['T2_C'] + d['T3_C']) / 2.0
        d['Tgap_abs'] = (d['T2_C'] - d['T3_C']).abs()
    else:
        # attempt to compute Tmean if one of T2 or T3 missing remains nan
        d['Tgap_abs'] = np.nan
    return d

def compute_interval_deltas(df: pd.DataFrame) -> pd.DataFrame:
    """
    For consecutive observations compute:
    - dt_s (interval duration)
    - dT2, dT3, dTmean
    - cooling rates (°C/h) for each variable (negative sign included as described)
    - interval_energy_removed_J and cumulative_energy_removed_J require mass and cp (handled elsewhere)
    """
    d = df.copy()
    d['dt_s'] = np.nan
    d['dT2'] = np.nan
    d['dT3'] = np.nan
    d['dTmean'] = np.nan
    d['rate_T2_C_per_h'] = np.nan
    d['rate_T3_C_per_h'] = np.nan
    d['rate_Tmean_C_per_h'] = np.nan

    if 'elapsed_s' not in d.columns:
        return d

    d.loc[:, 'dt_s'] = d['elapsed_s'].diff()
    if 'T2_C' in d.columns:
        d['dT2'] = d['T2_C'].diff()
        d['rate_T2_C_per_h'] = - d['dT2'] / d['dt_s'] * 3600.0
    if 'T3_C' in d.columns:
        d['dT3'] = d['T3_C'].diff()
        d['rate_T3_C_per_h'] = - d['dT3'] / d['dt_s'] * 3600.0
    if 'Tmean' in d.columns:
        d['dTmean'] = d['Tmean'].diff()
        d['rate_Tmean_C_per_h'] = - d['dTmean'] / d['dt_s'] * 3600.0

    # Handle zero or invalid dt_s by setting rates to NaN
    invalid = (d['dt_s'] <= 0) | (~np.isfinite(d['dt_s']))
    for c in ['rate_T2_C_per_h','rate_T3_C_per_h','rate_Tmean_C_per_h']:
        if c in d.columns:
            d.loc[invalid, c] = np.nan

    return d

def compute_interval_energy_and_cumulative(df: pd.DataFrame, mass_kg: float, cp_j_per_kgk: float) -> pd.DataFrame:
    """
    Interval energy removed (signed): m * cp * (Tmean_prev - Tmean_current)
    Cumulative sensible energy removed relative to first reading: m * cp * (Tmean_initial - Tmean_current)

    The first row will have interval energy NaN and cumulative energy zero.
    """
    d = df.copy()
    d['interval_energy_J'] = np.nan
    d['cumulative_energy_J'] = np.nan
    if 'Tmean' not in d.columns:
        return d
    m = float(mass_kg)
    cp = float(cp_j_per_kgk)
    if len(d) == 0:
        return d
    T0 = d['Tmean'].iloc[0]
    d['cumulative_energy_J'] = m * cp * (T0 - d['Tmean'])
    d['interval_energy_J'] = np.nan
    if len(d) >= 2:
        d['interval_energy_J'].iloc[1:] = m * cp * (d['Tmean'].iloc[:-1].values - d['Tmean'].iloc[1:].values)
    # Ensure dtype floats
    d['interval_energy_J'] = d['interval_energy_J'].astype(float)
    d['cumulative_energy_J'] = d['cumulative_energy_J'].astype(float)
    return d

def compute_interval_heat_loss_rate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Interval average net heat-loss rate (W) = interval_energy_J / dt_s
    Also compute average_net_heat_loss_rate_overall = total_cumulative_energy / total_elapsed_duration (caller)
    """
    d = df.copy()
    d['interval_heat_loss_W'] = np.nan
    if 'interval_energy_J' in d.columns and 'dt_s' in d.columns:
        d['interval_heat_loss_W'] = d['interval_energy_J'] / d['dt_s']
        # First row: NaN (interval_energy is NaN)
    return d

def average_net_heat_loss_rate_over_window(df: pd.DataFrame) -> Optional[float]:
    """
    Average net heat-loss rate (W) = total sensible energy removed (J) / total elapsed duration (s)
    Use cumulative energy at last reading and elapsed time difference between first and last valid readings.
    """
    if 'elapsed_s' not in df.columns or 'cumulative_energy_J' not in df.columns:
        return None
    if len(df) < 2:
        return None
    total_energy_removed = df['cumulative_energy_J'].iloc[-1]
    total_time = df['elapsed_s'].iloc[-1] - df['elapsed_s'].iloc[0]
    if total_time <= 0:
        return None
    return float(total_energy_removed / total_time)

def compute_apparent_heat_flux(avg_W: float, area_m2: float) -> Optional[float]:
    """
    Apparent average heat flux (W/m²) = avg_W / area_m2
    """
    if area_m2 is None or area_m2 == 0:
        return None
    return float(avg_W / area_m2)

# -----------------------
# Ambient-corrected
# -----------------------

def ambient_corrected_theta(df: pd.DataFrame, ambient_col: Optional[str]='T0_C') -> pd.DataFrame:
    """
    theta = Tmean - T0 (per-row); normalized theta = theta_current / theta_initial
    If ambient not available per-row, ambient_col may be None; caller should pass constant ambient separately.
    """
    d = df.copy()
    d['theta_C'] = np.nan
    d['theta_norm'] = np.nan
    if 'Tmean' not in d.columns:
        return d
    if ambient_col is not None and ambient_col in d.columns:
        d['theta_C'] = d['Tmean'] - d[ambient_col]
    else:
        # caller should set a constant ambient as another column, e.g., 'ambient_used_C'
        if 'ambient_used_C' in d.columns:
            d['theta_C'] = d['Tmean'] - d['ambient_used_C']
        else:
            # cannot compute
            return d
    # Normalized
    initial = d['theta_C'].iloc[0]
    if pd.isna(initial) or initial == 0:
        d['theta_norm'] = np.nan
    else:
        d['theta_norm'] = d['theta_C'] / initial
    return d

# -----------------------
# Fits & statistics
# -----------------------

def linear_fit_time(temp: np.ndarray, time_s: np.ndarray) -> Dict:
    """
    Fit Tmean = intercept + slope * t
    Returns dict with intercept, slope (°C/s), slope_per_min (°C/min), slope_per_hour (°C/h), R2, RMSE, y_fit, residuals
    """
    res = {}
    # drop nan pairs
    mask = (~np.isnan(temp)) & (~np.isnan(time_s))
    if mask.sum() < 2:
        return {'success': False, 'message': 'Insufficient points for linear fit'}
    x = time_s[mask].astype(float)
    y = temp[mask].astype(float)
    # linear polyfit
    p = np.polyfit(x, y, 1)
    slope = float(p[0])
    intercept = float(p[1])
    y_fit = intercept + slope * x
    residuals = y - y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan
    rmse = np.sqrt(np.mean(residuals**2))
    res.update({
        'success': True,
        'intercept': intercept,
        'slope_per_s': slope,
        'slope_per_min': slope * 60.0,
        'slope_per_h': slope * 3600.0,
        'R2': float(r2),
        'RMSE_C': float(rmse),
        'y_fit': (x, y_fit),
        'residuals': (x, residuals)
    })
    return res

def _exp_model(t_min, A, k):
    return A * np.exp(-k * t_min)

def exponential_fit_theta(theta: np.ndarray, time_s: np.ndarray) -> Dict:
    """
    Fit theta(t) = A * exp(-k * t) with t in minutes

    Returns dict: A, k (min^-1), tau (min), R2, RMSE, y_fit, residuals
    """
    res = {}
    mask = (~np.isnan(theta)) & (~np.isnan(time_s))
    if mask.sum() < 3:
        return {'success': False, 'message': 'Insufficient points for exponential fit'}
    tmin = time_s[mask].astype(float) / 60.0
    y = theta[mask].astype(float)
    # If y contains negative or zero values we'll still attempt to fit, but that usually fails — user must handle.
    # Provide initial guesses: A = y0, k = 0.01
    A0 = float(y[0]) if len(y) > 0 else 1.0
    k0 = 0.01
    try:
        popt, pcov = curve_fit(_exp_model, tmin, y, p0=[A0, k0], bounds=[[-np.inf, 0.0], [np.inf, np.inf]], maxfev=10000)
    except Exception as e:
        return {'success': False, 'message': f'Exponential fit failed: {e}'}
    A, k = float(popt[0]), float(popt[1])
    y_fit_vals = _exp_model(tmin, A, k)
    residuals = y - y_fit_vals
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan
    rmse = float(np.sqrt(np.mean(residuals**2)))
    tau = 1.0 / k if k > 0 else np.nan
    res.update({
        'success': True,
        'A': A,
        'k_min_inv': k,
        'tau_min': tau,
        'R2': float(r2),
        'RMSE_C': rmse,
        'y_fit': (tmin, y_fit_vals),
        'residuals': (tmin, residuals)
    })
    return res

# -----------------------
# Helper: threshold crossing
# -----------------------

def find_threshold_crossing(time_s: np.ndarray, value: np.ndarray, threshold: float, direction='down') -> Optional[float]:
    """
    Locate the first crossing of 'threshold' in value vs time. Only searches actual measurements and linearly interpolates
    between consecutive measurements if threshold lies between them. Does NOT extrapolate beyond recorded data.

    direction: 'down' (value decreases to cross threshold) or 'up' (value increases)
    Returns crossing time in seconds or None if not reached. If first reading already below threshold for 'down', returns 0.0
    """
    ts = np.asarray(time_s)
    vs = np.asarray(value)
    # require at least one non-nan
    valid = (~np.isnan(ts)) & (~np.isnan(vs))
    if valid.sum() == 0:
        return None
    ts = ts[valid]
    vs = vs[valid]
    # if first reading already meets criteria
    if direction == 'down':
        if vs[0] <= threshold:
            return float(0.0)
        # search for first pair where vs[i] >= threshold and vs[i+1] <= threshold
        for i in range(len(vs)-1):
            v1, v2 = vs[i], vs[i+1]
            if (v1 >= threshold and v2 <= threshold) or (v1 > threshold and v2 < threshold):
                # linear interpolation
                t1, t2 = ts[i], ts[i+1]
                if v1 == v2:
                    return float(t1)
                frac = (v1 - threshold) / (v1 - v2)
                return float(t1 + frac * (t2 - t1))
        return None
    elif direction == 'up':
        if vs[0] >= threshold:
            return float(0.0)
        for i in range(len(vs)-1):
            v1, v2 = vs[i], vs[i+1]
            if (v1 <= threshold and v2 >= threshold):
                t1, t2 = ts[i], ts[i+1]
                if v1 == v2:
                    return float(t1)
                frac = (threshold - v1) / (v2 - v1)
                return float(t1 + frac * (t2 - t1))
        return None
    else:
        raise ValueError("direction must be 'down' or 'up'")
