from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .constants import MOD_TO_BASE, VALID_BASES, base_for_label, canonical_base, normalize_sequence


@dataclass(frozen=True)
class ModelSpec:
    model_name: str
    positions: list[str]
    numeric: list[str]


BINARY_POSITION_SETS = {
    "m1": ["position1", "position2", "position3", "position4"],
    "m2": ["position2", "position3", "position4", "position6"],
    "m3": ["position3", "position4", "position6", "position7"],
    "m4": ["position4", "position6", "position7", "position8"],
    "m5": ["position6", "position7", "position8", "position9"],
}

MULTI_POSITION_SETS = {
    "m1": ["position1", "position2", "position3", "position4", "position5_ref"],
    "m2": ["position2", "position3", "position4", "position6", "position5_ref"],
    "m3": ["position3", "position4", "position6", "position7", "position5_ref"],
    "m4": ["position4", "position6", "position7", "position8", "position5_ref"],
    "m5": ["position6", "position7", "position8", "position9", "position5_ref"],
}


def metric_columns(metric: str) -> list[str]:
    return [f"{metric}_output{i}" for i in range(1, 6)]


def all_numeric_columns(condition: str) -> list[str]:
    if condition == "full":
        cols: list[str] = []
        for idx in range(1, 6):
            cols.extend(
                [
                    f"mean_output{idx}",
                    f"median_output{idx}",
                    f"std_output{idx}",
                    f"length_output{idx}",
                ]
            )
        return cols
    return metric_columns(condition)


def numeric_for_component(model_name: str, condition: str) -> list[str]:
    idx = int(model_name[1])
    if condition == "full":
        return [
            f"mean_output{idx}",
            f"median_output{idx}",
            f"std_output{idx}",
            f"length_output{idx}",
        ]
    return [f"{condition}_output{idx}"]


def binary_model_specs(strategy: str, condition: str) -> list[ModelSpec]:
    if strategy == "strategy1":
        return [
            ModelSpec(name, positions, numeric_for_component(name, condition))
            for name, positions in BINARY_POSITION_SETS.items()
        ]
    if strategy == "strategy2":
        nums = all_numeric_columns(condition)
        return [ModelSpec(name, positions, nums) for name, positions in BINARY_POSITION_SETS.items()]
    if strategy == "strategy3":
        return [
            ModelSpec(
                "m_all",
                [
                    "position1",
                    "position2",
                    "position3",
                    "position4",
                    "position6",
                    "position7",
                    "position8",
                    "position9",
                ],
                all_numeric_columns(condition),
            )
        ]
    raise ValueError(f"Unsupported strategy: {strategy}")


def multiclass_model_specs(strategy: str) -> list[ModelSpec]:
    if strategy == "strategy1":
        return [
            ModelSpec(name, positions, numeric_for_component(name, "full"))
            for name, positions in MULTI_POSITION_SETS.items()
        ]
    if strategy == "strategy2":
        nums = all_numeric_columns("full")
        return [ModelSpec(name, positions, nums) for name, positions in MULTI_POSITION_SETS.items()]
    if strategy == "strategy3":
        return [
            ModelSpec(
                "m_all",
                [
                    "position1",
                    "position2",
                    "position3",
                    "position4",
                    "position6",
                    "position7",
                    "position8",
                    "position9",
                    "position5_ref",
                ],
                all_numeric_columns("full"),
            )
        ]
    raise ValueError(f"Unsupported strategy: {strategy}")


def add_five_mer_positions(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence_input"] = df["sequence_input"].map(normalize_sequence)
    bad = df["sequence_input"].str.len() != 5
    if bad.any():
        raise ValueError(f"sequence_input must be 5 bases for all rows; bad rows: {bad.sum()}")
    for idx in range(5):
        df[f"position{idx + 1}"] = df["sequence_input"].str[idx]
    return df


def finalize_25mer(parts: list[str]) -> str | None:
    if len(parts) != 25:
        return None
    checks = (
        parts[1] == parts[5]
        and parts[2] == parts[6]
        and parts[3] == parts[7]
        and parts[4] == parts[8]
        and parts[6] == parts[10]
        and parts[7] == parts[11]
        and parts[8] == parts[12]
        and parts[9] == parts[13]
        and parts[11] == parts[15]
        and parts[12] == parts[16]
        and parts[13] == parts[17]
        and parts[14] == parts[18]
        and parts[16] == parts[20]
        and parts[17] == parts[21]
        and parts[18] == parts[22]
        and parts[19] == parts[23]
    )
    if not checks:
        return None
    return "_".join(parts[0:5] + parts[21:25])


def generate_9mer_candidates(extracted: pd.DataFrame) -> pd.DataFrame:
    """Generate centered 9-mer records from extracted 5-mer rows without sampling."""
    required = ["fast5_prefix", "fast5_order", "sequence_input"]
    missing = [c for c in required if c not in extracted.columns]
    if missing:
        raise ValueError(f"Cannot generate 9-mer candidates; missing columns: {missing}")

    df = add_five_mer_positions(extracted)
    df["fast5_order"] = pd.to_numeric(df["fast5_order"], errors="raise").astype(int)
    records: list[pd.DataFrame] = []

    for _, group in df.groupby("fast5_prefix", sort=False):
        group = group.sort_values("fast5_order").reset_index(drop=True)
        orders = group["fast5_order"].to_numpy()
        for start in range(len(group) - 4):
            window_orders = orders[start : start + 5]
            if not np.all(np.diff(window_orders) == 1):
                continue
            window = group.iloc[start : start + 5]
            parts: list[str] = []
            for _, row in window.iterrows():
                parts.extend([canonical_base(row[f"position{i}"]) for i in range(1, 6)])
            nine = "_".join(parts)
            final = finalize_25mer(parts)
            if final is None:
                continue
            center = parts[12]
            mid = window.iloc[[2]].copy()
            mid["9mer"] = nine
            mid["center"] = center
            mid["9mer_final"] = final
            mid["center_match"] = final.split("_")[4] == center
            mid["window_start_order"] = int(window_orders[0])
            mid["window_end_order"] = int(window_orders[-1])
            records.append(mid)

    if not records:
        raise ValueError("No valid 9-mer candidates were generated from the input.")
    result = pd.concat(records, ignore_index=True)
    return add_nine_mer_positions(result)


def add_nine_mer_positions(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "9mer_final" not in df.columns:
        raise ValueError("Prepared data must contain 9mer_final")
    parts = df["9mer_final"].astype(str).str.split("_", expand=True)
    if parts.shape[1] != 9:
        raise ValueError("9mer_final must split into exactly 9 positions")
    df["center_true_label"] = parts[4].astype(str).str.strip()
    for idx in range(9):
        df[f"position{idx + 1}"] = parts[idx].map(canonical_base)
    df["position5_ref"] = df["position5"].map(canonical_base)
    df["center_ref_base"] = df["position5_ref"]
    return df


def filter_binary_candidates(df: pd.DataFrame, base: str) -> pd.DataFrame:
    if "center_ref_base" not in df.columns:
        df = add_nine_mer_positions(df)
    return df.loc[df["center_ref_base"] == base].reset_index(drop=True)


def filter_by_target_bases(df: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    bases = {base_for_label(label) for label in labels}
    if "center_ref_base" not in df.columns:
        df = add_nine_mer_positions(df)
    return df.loc[df["center_ref_base"].isin(bases)].reset_index(drop=True)


def one_hot_positions(df: pd.DataFrame, pos_columns: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for col in pos_columns:
        canonical = df[col].map(canonical_base)
        for base in VALID_BASES:
            out[f"{col}_{base}"] = (canonical == base).astype("uint8")
    return out.reset_index(drop=True)


def build_feature_matrix(df: pd.DataFrame, spec: ModelSpec) -> pd.DataFrame:
    missing = [c for c in spec.positions + spec.numeric if c not in df.columns]
    if missing:
        raise ValueError(f"Missing feature columns for {spec.model_name}: {missing}")
    pos = one_hot_positions(df, spec.positions)
    nums = df[spec.numeric].apply(pd.to_numeric, errors="raise").reset_index(drop=True)
    return pd.concat([pos, nums], axis=1)


def validate_targets(labels: list[str], valid_labels: list[str]) -> list[str]:
    unknown = [label for label in labels if label not in valid_labels]
    if unknown:
        raise ValueError(f"Unknown target label(s): {unknown}. Valid labels: {valid_labels}")
    return labels
