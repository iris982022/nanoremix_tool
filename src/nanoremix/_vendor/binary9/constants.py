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
    "Am",
    "m5C",
    "hm5C",
    "ac4C",
    "o8G",
    "Gm",
    "m5U",
    "pseudoU",
    "Um",
    "m1psi",
]

TOKEN_ALIASES = {
    "T": "U",
    "normalA": "A",
    "normalC": "C",
    "normalG": "G",
    "normalU": "U",
    "I": "Inosine",
    "Ψ": "pseudoU",
    "psi": "pseudoU",
    "Psi": "pseudoU",
    "PSI": "pseudoU",
    "8-oxoG": "o8G",
    "8oxoG": "o8G",
    "oxoG": "o8G",
    "OxoG": "o8G",
    "OG": "o8G",
}

MOD_TO_BASE = {
    "A": "A",
    "m6A": "A",
    "m1A": "A",
    "Am": "A",
    "Inosine": "A",
    "C": "C",
    "m5C": "C",
    "hm5C": "C",
    "FormylC": "C",
    "ac4C": "C",
    "Cm": "C",
    "G": "G",
    "m7G": "G",
    "o8G": "G",
    "Gm": "G",
    "U": "U",
    "pseudoU": "U",
    "m1psi": "U",
    "m5U": "U",
    "Um": "U",
}

DEFAULT_MODEL_ROOT = Path("models/RNA004/binary")

DEFAULT_MULTICLASS_MODEL_ROOT = Path("models/RNA004/multi")

MULTICLASS_STRATEGY_DIR_NAMES = {
    "strategy1": "ref_4current_8pos",
    "strategy2": "ref_all_8pos",
    "strategy3": "ref_all_16pos",
}

DEFAULT_MULTICLASS_MODEL_DIRS = {
    strategy: DEFAULT_MULTICLASS_MODEL_ROOT / dirname
    for strategy, dirname in MULTICLASS_STRATEGY_DIR_NAMES.items()
}

MULTICLASS_STRATEGY_ENV_VARS = {
    "strategy1": "NANOREMIX9MER_MULTICLASS_STRATEGY1_DIR",
    "strategy2": "NANOREMIX9MER_MULTICLASS_STRATEGY2_DIR",
    "strategy3": "NANOREMIX9MER_MULTICLASS_STRATEGY3_DIR",
}


def normalize_token(token: object) -> str:
    value = str(token).strip()
    return TOKEN_ALIASES.get(value, value)


def canonical_base(token: object) -> str:
    value = normalize_token(token)
    if value not in MOD_TO_BASE:
        raise ValueError(f"Unknown base/modification token: {value}")
    return MOD_TO_BASE[value]


def normalize_sequence(seq: object) -> str:
    return str(seq).upper().replace("T", "U")


def base_for_label(label: str) -> str:
    value = normalize_token(label)
    if value not in MOD_TO_BASE:
        raise ValueError(f"Unknown modification label: {label}")
    return MOD_TO_BASE[value]
