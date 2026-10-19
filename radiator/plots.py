"""
Plotly plots used by app.py
"""

import plotly.graph_objs as go
import pandas as pd
import numpy as np

def plot_temperature(df: pd.DataFrame, title="Temperature vs Time", fit_curve=None):
    fig = go.Figure()
    t_min = df['elapsed_s'] / 60.0
    if 'T0_C' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['T0_C'], mode='lines+markers', name='T0 (ambient °C)'))
    if 'T1_C' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['T1_C'], mode='lines+markers', name='T1 (heated water °C)'))
    if 'T2_C' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['T2_C'], mode='lines+markers', name='T2 inlet °C'))
    if 'T3_C' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['T3_C'], mode='lines+markers', name='T3 outlet °C'))
    if 'Tmean' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['Tmean'], mode='lines+markers', name='Tmean (endpoints mean °C)'))
    if fit_curve is not None and fit_curve.get('x') is not None:
        fig.add_trace(go.Scatter(x=np.array(fit_curve['x'])/60.0, y=fit_curve['y'], mode='lines', name='Fitted curve', line=dict(dash='dash')))
    fig.update_layout(title=title, xaxis_title='Elapsed time (min)', yaxis_title='Temperature (°C)', hovermode='closest', legend=dict(orientation='h'))
    return fig

def plot_theta(df: pd.DataFrame, title="Ambient-corrected temperature (θ)"):
    fig = go.Figure()
    if 'theta_C' in df.columns:
        t_min = df['elapsed_s'] / 60.0
        fig.add_trace(go.Scatter(x=t_min, y=df['theta_C'], mode='lines+markers', name='θ = Tmean - T0 (°C)'))
    fig.update_layout(title=title, xaxis_title='Elapsed time (min)', yaxis_title='θ (°C)', hovermode='closest', legend=dict(orientation='h'))
    return fig

def plot_energy_cumulative(df: pd.DataFrame, title="Cumulative sensible energy removed"):
    fig = go.Figure()
    if 'cumulative_energy_J' in df.columns:
        t_min = df['elapsed_s'] / 60.0
        fig.add_trace(go.Scatter(x=t_min, y=df['cumulative_energy_J'], mode='lines+markers', name='Cumulative energy (J)'))
    fig.update_layout(title=title, xaxis_title='Elapsed time (min)', yaxis_title='Energy (J)', hovermode='closest', legend=dict(orientation='h'))
    return fig

def plot_cooling_rates(df: pd.DataFrame, title="Interval cooling rates"):
    fig = go.Figure()
    t_min = df['elapsed_s'] / 60.0
    if 'rate_T2_C_per_h' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['rate_T2_C_per_h'], mode='lines+markers', name='Rate T2 (°C/h)'))
    if 'rate_T3_C_per_h' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['rate_T3_C_per_h'], mode='lines+markers', name='Rate T3 (°C/h)'))
    if 'rate_Tmean_C_per_h' in df.columns:
        fig.add_trace(go.Scatter(x=t_min, y=df['rate_Tmean_C_per_h'], mode='lines+markers', name='Rate Tmean (°C/h)'))
    fig.update_layout(title=title, xaxis_title='Elapsed time (min)', yaxis_title='Cooling rate (°C/h)', hovermode='closest', legend=dict(orientation='h'))
    return fig
