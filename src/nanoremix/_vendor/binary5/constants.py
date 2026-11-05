from __future__ import annotations

from pathlib import Path

VALID_BASES = ("A", "C", "G", "U")
VALID_STRATEGIES = ("strategy1", "strategy2", "strategy3")
VALID_BINARY_CONDITIONS = ("full", "mean", "std", "length")

BINARY_PAIRS = {
    "A_m6A": {"base": "A", "mod": "m6A"},
    "C_m5C": {"base": "C", "mod": "m5C"},
    "U_pseudoU": {"base": "U", "mod": "pseudoU"},
}

MULTICLASS_LABELS = [
    "A",
    "C",
    "G",
    "U",
    "m6A",
    "m1A",
    "Inosine",
    "Am",
    "m5C",
    "hm5C",
    "ac4C",
    "Cm",
    "FormylC",
    "pseudoU",
    "m1psi",
    "m5U",
    "Um",
    "m7G",
    "o8G",
]

MOD_TO_BASE = {
    "A": "A",
    "normalA": "A",
    "m6A": "A",
    "m1A": "A",
    "Inosine": "A",
    "Am": "A",
    "C": "C",
    "normalC": "C",
    "m5C": "C",
    "hm5C": "C",
    "FormylC": "C",
    "ac4C": "C",
    "Cm": "C",
    "G": "G",
    "normalG": "G",
    "m7G": "G",
    "o8G": "G",
    "Gm": "G",
    "U": "U",
    "T": "U",
    "normalU": "U",
    "pseudoU": "U",
    "m1psi": "U",
    "m5U": "U",
    "Um": "U",
}

DEFAULT_MODEL_ROOT = Path("models/RNA002/binary")


def canonical_base(token: object) -> str:
    value = str(token)
    if value not in MOD_TO_BASE:
        raise ValueError(f"Unknown base/modification token: {value}")
    return MOD_TO_BASE[value]


def normalize_sequence(seq: object) -> str:
    return str(seq).upper().replace("T", "U")


def base_for_label(label: str) -> str:
    if label not in MOD_TO_BASE:
        raise ValueError(f"Unknown modification label: {label}")
    return MOD_TO_BASE[label]
