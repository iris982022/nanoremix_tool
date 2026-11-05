from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import normalize_sequence

METRIC_PREFIXES = ("mean", "std", "median", "length")


def numeric_columns() -> list[str]:
    cols: list[str] = []
    for metric in METRIC_PREFIXES:
        cols.extend([f"{metric}_output{i}" for i in range(1, 10)])
    return cols


STEP6_REQUIRED_COLUMNS = ["read_id", "ref_start", "sequence_input", *numeric_columns()]


def _table_separator(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def coerce_feature_types(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "sequence_input" in df.columns:
        df["sequence_input"] = df["sequence_input"].map(normalize_sequence)
        bad = df["sequence_input"].str.len() != 9
        if bad.any():
            raise ValueError(f"sequence_input must be 9 bases for all rows; bad rows: {int(bad.sum())}")
    if "ref_start" in df.columns:
        df["ref_start"] = pd.to_numeric(df["ref_start"], errors="raise").astype(int)
    for col in numeric_columns():
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="raise")
    return df


def read_step6_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path, sep=_table_separator(path), low_memory=False)
    missing = [c for c in STEP6_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"step6 input is missing required columns: {missing}")
    return coerce_feature_types(df)


def read_table_auto(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    return pd.read_csv(path, sep=_table_separator(path), low_memory=False)


def detect_input_format(path: str | Path) -> str:
    path = Path(path)
    header = pd.read_csv(path, sep=None, engine="python", nrows=0)
    cols = set(header.columns)
    if "9mer_final" in cols:
        return "prepared"
    if {"read_id", "ref_start", "sequence_input"}.issubset(cols):
        return "step6"
    raise ValueError("Could not auto-detect input format. Use --input-format step6 or prepared.")
