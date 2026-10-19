"""
I/O helpers: CSV loading, column alias detection and validation, processed-data exports.
"""

import pandas as pd
from typing import Tuple, List, Dict
from radiator.calculations import map_columns, prepare_dataframe

REQUIRED_TEMPS = ['T2_C','T3_C']

def load_csv_file(path_or_buffer) -> Tuple[pd.DataFrame, dict]:
    """
    Load a CSV into a DataFrame. Return df and detection mapping
    """
    df = pd.read_csv(path_or_buffer)
    mapping = map_columns(df)
    return df, mapping

def validate_and_prepare(df: pd.DataFrame, mapping: dict, require_elapsed=True) -> Tuple[pd.DataFrame, List[str]]:
    """
    Returns processed dataframe ready for calculations and a list of messages describing missing/available columns.
    """
    messages = []
    mapped = mapping.copy()
    # Prepare dataframe including creation of elapsed_s if datetime present
    d = prepare_dataframe(df, mapped, require_elapsed=require_elapsed)
    # Report which required columns are present
    for r in REQUIRED_TEMPS:
        if r not in d.columns or d[r].isna().all():
            messages.append(f"Missing or empty required temperature column: {r}")
        else:
            messages.append(f"Found {r}")
    if 'elapsed_s' not in d.columns or d['elapsed_s'].isna().all():
        messages.append("Missing elapsed_s column or invalid; elapsed time required for analysis.")
    else:
        messages.append("Found elapsed_s")
    return d, messages
