import math
import numpy as np
import pandas as pd
from radiator import calculations as calc

def test_volume_to_mass():
    mass = calc.volume_ml_to_mass_kg(1000.0, density_kg_m3=1000.0)  # 1 L of water at density 1000
    assert pytest_approx(mass, 1.0, rel=1e-6)

def test_rectangular_area():
    area = calc.rectangular_tube_area_m2(length_mm=100.0, width_mm=10.0, height_mm=5.0)
    # area = 2*(0.01+0.005)*0.1 = 2*(0.015)*0.1 = 0.003 m2
    assert pytest_approx(area, 0.003, rel=1e-6)

def test_energy_and_heat_loss():
    # synthetic data: Tmean starts 60 and goes to 50 in 10s, mass 0.02 kg cp 4180
    df = pd.DataFrame({'elapsed_s':[0,10], 'T2_C':[61,51], 'T3_C':[59,49]})
    df = calc.compute_tmean_and_gap(df)
    df = calc.compute_interval_deltas(df)
    mass = calc.volume_ml_to_mass_kg(20.0, 998.0)
    df = calc.compute_interval_energy_and_cumulative(df, mass, 4180.0)
    df = calc.compute_interval_heat_loss_rate(df)
    # cumulative energy = m*cp*(T0 - Tcurrent): initial - final = m*cp*(60 - 50) = m*cp*10
    expected = mass * 4180.0 * ((df['Tmean'].iloc[0] - df['Tmean'].iloc[-1]))
    assert pytest_approx(df['cumulative_energy_J'].iloc[-1], expected, rel=1e-6)
    # interval heat loss rate = interval_energy / dt = expected / 10
    assert pytest_approx(df['interval_heat_loss_W'].iloc[-1], expected / 10.0, rel=1e-6)

def test_cooling_rate_signs():
    # decreasing temperature yields positive cooling rate (°C/h)
    df = pd.DataFrame({'elapsed_s':[0,3600], 'T2_C':[60,59], 'T3_C':[60,59]})
    df = calc.compute_tmean_and_gap(df)
    df = calc.compute_interval_deltas(df)
    # dTmean = -1, dt = 3600 -> rate = -(-1)/3600*3600 = 1 °C/h
    assert pytest_approx(df['rate_Tmean_C_per_h'].iloc[-1], 1.0, rel=1e-6)

def test_linear_and_exponential_fit():
    # Linear: T = 100 - 0.5 * t_min ; use seconds
    t_s = np.array([0,60,120,180], dtype=float)
    t_min = t_s / 60.0
    T = 100.0 - 0.5 * t_min
    lin = calc.linear_fit_time(T, t_s)
    assert lin['success']
    assert pytest_approx(lin['slope_per_min'], -0.5, rel=1e-3)
    # Exponential: theta = A*exp(-k t_min)
    A = 20.0
    k = 0.1
    theta = A * np.exp(-k * t_min)
    expres = calc.exponential_fit_theta(theta, t_s)
    assert expres['success']
    assert pytest_approx(expres['A'], A, rel=1e-2)
    assert pytest_approx(expres['k_min_inv'], k, rel=1e-2)

def test_threshold_crossing():
    t = np.array([0,10,20,30])
    val = np.array([100,80,60,40])
    # threshold 70 should be crossed between 10 and 20 -> interpolation: 10 + (80-70)/(80-60)*10 = 10 + 0.5*10 = 15
    tcross = calc.find_threshold_crossing(t, val, 70.0, direction='down')
    assert pytest_approx(tcross, 15.0, rel=1e-6)

# helper approx using pytest style without importing pytest floats
def pytest_approx(a,b,rel=1e-6):
    return abs(a-b) <= max(rel*abs(b), 1e-9)
