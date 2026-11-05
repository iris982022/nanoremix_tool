from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .constants import VALID_BASES, base_for_label, canonical_base, normalize_sequence, normalize_token


@dataclass(frozen=True)
class ModelSpec:
    model_name: str
    positions: list[str]
    numeric: list[str]


COMPONENT_POSITION_SETS = {
    "m1": [
        "position1",
        "position2",
        "position3",
        "position4",
        "position5",
        "position6",
        "position7",
        "position8",
    ],
    "m2": [
        "position2",
        "position3",
        "position4",
        "position5",
        "position6",
        "position7",
        "position8",
        "position10",
    ],
    "m3": [
        "position3",
        "position4",
        "position5",
        "position6",
        "position7",
        "position8",
        "position10",
        "position11",
    ],
    "m4": [
        "position4",
        "position5",
        "position6",
        "position7",
        "position8",
        "position10",
        "position11",
        "position12",
    ],
    "m5": [
        "position5",
        "position6",
        "position7",
        "position8",
        "position10",
        "position11",
        "position12",
        "position13",
    ],
    "m6": [
        "position6",
        "position7",
        "position8",
        "position10",
        "position11",
        "position12",
        "position13",
        "position14",
    ],
    "m7": [
        "position7",
        "position8",
        "position10",
        "position11",
        "position12",
        "position13",
        "position14",
        "position15",
    ],
    "m8": [
        "position8",
        "position10",
        "position11",
        "position12",
        "position13",
        "position14",
        "position15",
        "position16",
    ],
    "m9": [
        "position10",
        "position11",
        "position12",
        "position13",
        "position14",
        "position15",
        "position16",
        "position17",
    ],
}

ALL_FLANK_POSITIONS = [
    "position1",
    "position2",
    "position3",
    "position4",
    "position5",
    "position6",
    "position7",
    "position8",
    "position10",
    "position11",
    "position12",
    "position13",
    "position14",
    "position15",
    "position16",
    "position17",
]


def metric_columns(metric: str) -> list[str]:
    return [f"{metric}_output{i}" for i in range(1, 10)]


def all_numeric_columns(condition: str) -> list[str]:
    if condition == "full":
        cols: list[str] = []
        for idx in range(1, 10):
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
    idx = int(model_name[1:])
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
            for name, positions in COMPONENT_POSITION_SETS.items()
        ]
    if strategy == "strategy2":
        nums = all_numeric_columns(condition)
        return [ModelSpec(name, positions, nums) for name, positions in COMPONENT_POSITION_SETS.items()]
    if strategy == "strategy3":
        return [ModelSpec("m_all", ALL_FLANK_POSITIONS, all_numeric_columns(condition))]
    raise ValueError(f"Unsupported strategy: {strategy}")


def multiclass_model_specs(strategy: str) -> list[ModelSpec]:
    if strategy == "strategy1":
        return [
            ModelSpec(name, positions + ["position9_ref"], numeric_for_component(name, "full"))
            for name, positions in COMPONENT_POSITION_SETS.items()
        ]
    if strategy == "strategy2":
        nums = all_numeric_columns("full")
        return [
            ModelSpec(name, positions + ["position9_ref"], nums)
            for name, positions in COMPONENT_POSITION_SETS.items()
        ]
    if strategy == "strategy3":
        return [ModelSpec("m_all", ALL_FLANK_POSITIONS + ["position9_ref"], all_numeric_columns("full"))]
    raise ValueError(f"Unsupported strategy: {strategy}")


def _sequence_tokens(seq: object) -> list[str]:
    value = normalize_sequence(seq)
    if len(value) != 9:
        raise ValueError(f"sequence_input must be 9 bases, observed {len(value)}: {value}")
    return [canonical_base(base) for base in value]


def _overlap_is_valid(rows: list[list[str]]) -> bool:
    return all(rows[idx][1:] == rows[idx + 1][:-1] for idx in range(len(rows) - 1))


def generate_17mer_candidates(step6: pd.DataFrame) -> pd.DataFrame:
    """Generate centered 9-mer model candidates from step6 data_df.csv rows."""
    required = ["read_id", "ref_start", "sequence_input"]
    missing = [c for c in required if c not in step6.columns]
    if missing:
        raise ValueError(f"Cannot generate 17-position candidates; missing columns: {missing}")

    df = step6.copy()
    df["sequence_input"] = df["sequence_input"].map(normalize_sequence)
    df["ref_start"] = pd.to_numeric(df["ref_start"], errors="raise").astype(int)
    df["_tokens"] = df["sequence_input"].map(_sequence_tokens)

    records: list[pd.DataFrame] = []
    for _, group in df.groupby("read_id", sort=False):
        group = group.sort_values("ref_start").reset_index(drop=True)
        starts = group["ref_start"].to_numpy(dtype=int)
        tokens = group["_tokens"].tolist()
        for start in range(len(group) - 8):
            window_starts = starts[start : start + 9]
            if not np.all(np.diff(window_starts) == 1):
                continue
            window_tokens = tokens[start : start + 9]
            if not _overlap_is_valid(window_tokens):
                continue

            tokens81 = [token for row_tokens in window_tokens for token in row_tokens]
            final17 = window_tokens[0] + window_tokens[-1][1:]
            center = window_tokens[4][4]

            mid = group.iloc[[start + 4]].drop(columns=["_tokens"]).copy()
            mid["9mer"] = "_".join(tokens81)
            mid["center"] = center
            mid["9mer_final"] = "_".join(final17)
            mid["center_match"] = center == final17[8]
            mid["window_start_ref_start"] = int(window_starts[0])
            mid["window_end_ref_start"] = int(window_starts[-1])
            records.append(mid.iloc[0].to_dict())

    if not records:
        raise ValueError("No valid 9-row windows were generated from the input.")
    result = pd.DataFrame.from_records(records)
    return add_seventeen_mer_positions(result)


def add_seventeen_mer_positions(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "9mer_final" not in df.columns:
        raise ValueError("Prepared data must contain 9mer_final")
    parts = df["9mer_final"].astype(str).str.split("_", expand=True)
    if parts.shape[1] != 17:
        raise ValueError("9mer_final must split into exactly 17 positions")
    df["center_true_label"] = parts[8].map(normalize_token)
    for idx in range(17):
        df[f"position{idx + 1}"] = parts[idx].map(normalize_token)
    df["position9_ref"] = df["position9"].map(canonical_base)
    df["center"] = df["position9"]
    df["center_ref_base"] = df["position9_ref"]
    return df


def filter_binary_candidates(df: pd.DataFrame, base: str) -> pd.DataFrame:
    if "center_ref_base" not in df.columns:
        df = add_seventeen_mer_positions(df)
    return df.loc[df["center_ref_base"] == base].reset_index(drop=True)


def filter_by_target_bases(df: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    bases = {base_for_label(label) for label in labels}
    if "center_ref_base" not in df.columns:
        df = add_seventeen_mer_positions(df)
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
