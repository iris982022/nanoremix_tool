from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from .constants import normalize_sequence

METRIC_PREFIXES = ("mean", "std", "median", "length")
EXTRACTED_COLUMNS = [
    "id",
    "sequence_input",
    "mean_output1",
    "mean_output2",
    "mean_output3",
    "mean_output4",
    "mean_output5",
    "std_output1",
    "std_output2",
    "std_output3",
    "std_output4",
    "std_output5",
    "median_output1",
    "median_output2",
    "median_output3",
    "median_output4",
    "median_output5",
    "length_output1",
    "length_output2",
    "length_output3",
    "length_output4",
    "length_output5",
    "fast5_prefix",
    "fast5_order",
]


def split_pipe_metric(value: object, metric_name: str) -> list[float]:
    parts = str(value).split("|")
    if len(parts) != 5:
        raise ValueError(f"{metric_name} must contain 5 pipe-separated values: {value}")
    return [float(x) for x in parts]


def parse_read_id(value: object) -> tuple[str, int | None, str | None]:
    parts = str(value).split("|")
    prefix = parts[0]
    order = None
    strand = None
    if len(parts) >= 2:
        try:
            order = int(parts[1])
        except ValueError:
            order = None
    if len(parts) >= 3:
        strand = parts[2]
    return prefix, order, strand


def read_step4_tsv(path: str | Path) -> pd.DataFrame:
    """Read NanoRemix step4 feature TSV and return extracted_data-style columns."""
    path = Path(path)
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, low_memory=False)
    if raw.shape[1] < 7:
        raise ValueError(
            f"Expected at least 7 columns in step4 TSV, observed {raw.shape[1]}: {path}"
        )

    first_id = str(raw.iloc[0, 0]).strip().lower()
    first_metric = str(raw.iloc[0, 3]).strip().lower() if raw.shape[1] > 3 else ""
    if first_id in {"id", "read_id", "readid"} or first_metric in {"mean", "mean_output"}:
        raw = raw.iloc[1:].reset_index(drop=True)

    raw = raw.iloc[:, :8].copy() if raw.shape[1] >= 8 else raw.iloc[:, :7].copy()
    cols = ["id", "base", "sequence_input", "mean", "std", "median", "length"]
    if raw.shape[1] == 8:
        cols.append("kmer_raw")
    raw.columns = cols

    records: dict[str, Iterable] = {
        "id": raw["id"].astype(str),
        "sequence_input": raw["sequence_input"].map(normalize_sequence),
    }

    for metric in METRIC_PREFIXES:
        split = raw[metric].map(lambda value, name=metric: split_pipe_metric(value, name))
        for idx in range(5):
            records[f"{metric}_output{idx + 1}"] = split.map(lambda values, i=idx: values[i])

    parsed_ids = raw["id"].map(parse_read_id)
    records["fast5_prefix"] = parsed_ids.map(lambda item: item[0])
    records["fast5_order"] = parsed_ids.map(lambda item: item[1])
    records["strand"] = parsed_ids.map(lambda item: item[2])
    if "base" in raw:
        records["raw_base"] = raw["base"].map(normalize_sequence)

    df = pd.DataFrame(records)
    if df["fast5_order"].isna().any():
        df["fast5_order"] = df.groupby("fast5_prefix").cumcount()
    df["fast5_order"] = df["fast5_order"].astype(int)
    return coerce_feature_types(df)


def coerce_feature_types(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence_input"] = df["sequence_input"].map(normalize_sequence)
    for metric in METRIC_PREFIXES:
        for idx in range(1, 6):
            col = f"{metric}_output{idx}"
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="raise")
    if "fast5_order" in df.columns:
        df["fast5_order"] = pd.to_numeric(df["fast5_order"], errors="raise").astype(int)
    return df


def read_extracted_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep, low_memory=False)
    missing = [c for c in EXTRACTED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Extracted input is missing required columns: {missing}")
    return coerce_feature_types(df)


def read_table_auto(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def detect_input_format(path: str | Path) -> str:
    path = Path(path)
    sample = pd.read_csv(path, sep=None, engine="python", header=None, nrows=2, dtype=str)
    if sample.shape[1] >= 7 and "|" in str(sample.iloc[0, 3]):
        return "step4"
    header = pd.read_csv(path, sep=None, engine="python", nrows=0)
    cols = set(header.columns)
    if "9mer_final" in cols:
        return "prepared"
    if {"sequence_input", "fast5_prefix", "fast5_order"}.issubset(cols):
        return "extracted"
    raise ValueError(
        "Could not auto-detect input format. Use --input-format step4, extracted, or prepared."
    )
