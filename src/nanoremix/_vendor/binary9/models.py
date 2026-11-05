from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import (
    DEFAULT_MODEL_ROOT,
    DEFAULT_MULTICLASS_MODEL_DIRS,
    MULTICLASS_LABELS,
    MULTICLASS_STRATEGY_DIR_NAMES,
    MULTICLASS_STRATEGY_ENV_VARS,
)


class LocalXGB:
    """Compatibility shim for 9mer multiclass pkls saved from __main__.LocalXGB."""

    _estimator_type = "classifier"

    def __init__(self, base_params=None, K=None, present_classes=None, num_global=None):
        self.base_params = dict(base_params or {})
        self.K = int(K or num_global or len(MULTICLASS_LABELS))
        present = [] if present_classes is None else present_classes
        self.present_ = np.array(sorted(int(c) for c in np.unique(present)))
        self.present = self.present_
        self.model_ = None
        self.model = None
        self.classes_ = np.arange(self.K)

    def fit(self, X, y, sample_weight=None, seed=0, verbose=False):
        from xgboost import XGBClassifier

        present = np.array(sorted(np.unique(y)))
        mapping = {c: i for i, c in enumerate(present)}
        y_local = pd.Series(y).map(mapping).astype(int).values
        params = dict(self.base_params)
        params.update(
            {
                "objective": "multi:softprob",
                "num_class": len(present),
                "eval_metric": "mlogloss",
                "tree_method": params.get("tree_method", "hist"),
                "random_state": params.get("random_state", seed),
            }
        )
        self.model_ = XGBClassifier(**params)
        self.model_.fit(X, y_local, sample_weight=sample_weight, verbose=verbose)
        self.present_ = present
        self.present = present
        return self

    def predict_proba(self, X):
        model = self.model_ if self.model_ is not None else self.model
        if model is None:
            raise ValueError("LocalXGB model is not loaded")
        local = model.predict_proba(X)
        if local.ndim == 1:
            local = np.vstack([1 - local, local]).T
        full = np.zeros((local.shape[0], self.K), dtype=float)
        present = self.present_ if self.present_ is not None and len(self.present_) else self.present
        for local_idx, global_idx in enumerate(present):
            if 0 <= int(global_idx) < self.K:
                full[:, int(global_idx)] = local[:, local_idx]
        denom = full.sum(axis=1, keepdims=True)
        denom[denom == 0] = 1.0
        return full / denom


def install_pickle_shims() -> None:
    main_mod = sys.modules.get("__main__")
    if main_mod is not None and not hasattr(main_mod, "LocalXGB"):
        setattr(main_mod, "LocalXGB", LocalXGB)



def _repair_xgboost_pickle_compat(root):
    """Fill attributes absent from XGBClassifier objects saved by older XGBoost."""
    try:
        from xgboost import XGBClassifier
        defaults = XGBClassifier().get_params(deep=False)
    except Exception:
        return

    defaults.setdefault("use_label_encoder", False)
    seen = set()

    def visit(obj):
        if obj is None or id(obj) in seen:
            return
        seen.add(id(obj))

        if isinstance(obj, XGBClassifier):
            for name, default in defaults.items():
                if not hasattr(obj, name):
                    setattr(obj, name, default)
            return

        if isinstance(obj, dict):
            for value in obj.values():
                visit(value)
            return

        if isinstance(obj, (list, tuple, set)):
            for value in obj:
                visit(value)
            return

        for attr in (
            "calibrated_classifiers_",
            "estimator",
            "base_estimator",
            "estimators_",
            "model_",
            "model",
            "classifier",
            "named_steps",
        ):
            if hasattr(obj, attr):
                visit(getattr(obj, attr))

    visit(root)

def load_model(path: str | Path):
    install_pickle_shims()
    try:
        import joblib
    except ImportError as exc:
        raise RuntimeError(
            "joblib is required to load NanoRemix model pkls. Install with `pip install -e .`."
        ) from exc
    model = joblib.load(path)
    _repair_xgboost_pickle_compat(model)
    return model


def _model_names_for_strategy(strategy: str) -> list[str]:
    return ["m_all"] if strategy == "strategy3" else [f"m{i}" for i in range(1, 10)]


def _paths_for_model_dir(model_dir: Path, strategy: str) -> dict[str, Path]:
    return {
        name: model_dir / f"{name}_calibrated.pkl"
        for name in _model_names_for_strategy(strategy)
    }


def _missing_model_paths(model_dir: Path, strategy: str) -> list[str]:
    return [
        str(path)
        for path in _paths_for_model_dir(model_dir, strategy).values()
        if not path.is_file()
    ]


def normalize_model_root(model_root: str | Path | None) -> Path:
    raw = model_root or os.environ.get("NANOREMIX9MER_MODEL_ROOT") or DEFAULT_MODEL_ROOT
    root = Path(raw).expanduser().resolve()
    if (root / "binary_classify").is_dir() and (
        (root / "multi_classify").is_dir() or (root / "mutli_classify").is_dir()
    ):
        return root
    if root.name in {"binary_classify", "multi_classify", "mutli_classify"}:
        return root.parent
    if (root / "9mer" / "binary_classify").is_dir():
        return root / "9mer"
    raise FileNotFoundError(
        f"Cannot locate NanoRemix 9mer model root from {root}. "
        "Pass --model-root to the directory containing binary_classify and multi_classify."
    )


def _find_strategy_dir(parent: Path, strategy: str) -> Path:
    candidate = parent / strategy
    if candidate.is_dir():
        return candidate
    for child in parent.iterdir():
        if child.is_dir() and child.name.strip() == strategy:
            return child
    raise FileNotFoundError(f"Strategy directory not found: {parent}/{strategy}")


def _multiclass_override_candidates(
    model_root: str | Path | None,
    strategy: str,
) -> list[Path]:
    candidates: list[Path] = []

    if model_root is not None:
        root = Path(model_root)
        mapped_dir = MULTICLASS_STRATEGY_DIR_NAMES[strategy]
        candidates.extend(
            [root / "multi_classify" / strategy, root / "mutli_classify" / strategy,
             root / mapped_dir, root / strategy]
        )
        if any(
            (root / f"{name}_calibrated.pkl").exists()
            for name in _model_names_for_strategy(strategy)
        ):
            candidates.append(root)

    env_var = MULTICLASS_STRATEGY_ENV_VARS[strategy]
    if os.environ.get(env_var):
        candidates.append(Path(os.environ[env_var]))

    mapped_dir = MULTICLASS_STRATEGY_DIR_NAMES[strategy]
    multi_root = os.environ.get("NANOREMIX9MER_MULTICLASS_MODEL_ROOT")
    if multi_root:
        root = Path(multi_root)
        candidates.extend([root / mapped_dir, root / strategy, root])

    candidates.append(DEFAULT_MULTICLASS_MODEL_DIRS[strategy])

    deduped: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        expanded = candidate.expanduser()
        key = str(expanded)
        if key not in seen:
            deduped.append(expanded)
            seen.add(key)
    return deduped


def _resolve_multiclass_model_dir(model_root: str | Path | None, strategy: str) -> Path:
    incomplete: list[tuple[Path, list[str]]] = []
    for candidate in _multiclass_override_candidates(model_root, strategy):
        if not candidate.is_dir():
            continue
        missing = _missing_model_paths(candidate, strategy)
        if not missing:
            return candidate
        incomplete.append((candidate, missing))

    if incomplete:
        details = "; ".join(
            f"{model_dir}: missing {missing}" for model_dir, missing in incomplete
        )
        raise FileNotFoundError(
            "Found configured NanoRemix 9mer multiclass model dir(s), "
            f"but the pkl set is incomplete: {details}"
        )

    root = normalize_model_root(model_root)
    multi_root = root / "multi_classify"
    if not multi_root.is_dir():
        multi_root = root / "mutli_classify"
    return _find_strategy_dir(multi_root, strategy)


def resolve_model_paths(
    model_root: str | Path | None,
    mode: str,
    strategy: str,
    condition: str = "full",
    pair: str | None = None,
) -> dict[str, Path]:
    if mode == "binary":
        root = normalize_model_root(model_root)
        if pair is None:
            raise ValueError("--pair is required for binary mode")
        base_dir = root / "binary_classify" / pair
        strategy_dir = _find_strategy_dir(base_dir, strategy)
        model_dir = strategy_dir if condition == "full" else strategy_dir / f"only_{condition}"
    elif mode == "multiclass":
        if condition != "full":
            raise ValueError("9mer multiclass pkls are fixed full-feature models; use --condition full.")
        model_dir = _resolve_multiclass_model_dir(model_root, strategy)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    paths = _paths_for_model_dir(model_dir, strategy)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing model pkl(s): {missing}")
    return paths


def feature_names_in(model) -> list[str] | None:
    names = getattr(model, "feature_names_in_", None)
    if names is not None:
        return list(names)
    calibrated = getattr(model, "calibrated_classifiers_", None)
    if calibrated:
        estimator = getattr(calibrated[0], "estimator", None) or getattr(
            calibrated[0], "base_estimator", None
        )
        names = getattr(estimator, "feature_names_in_", None)
        if names is not None:
            return list(names)
        inner = getattr(estimator, "model_", None) or getattr(estimator, "model", None)
        names = getattr(inner, "feature_names_in_", None)
        if names is not None:
            return list(names)
    return None


def align_feature_matrix(X, model):
    names = feature_names_in(model)
    if not names:
        return X
    missing = [name for name in names if name not in X.columns]
    extra = [name for name in X.columns if name not in names]
    if missing or extra:
        raise ValueError(
            "Feature columns do not match the trained model. "
            f"Missing: {missing}; extra: {extra}"
        )
    return X[names]


def binary_positive_probability(model, X) -> np.ndarray:
    prob = model.predict_proba(X)
    classes = list(getattr(model, "classes_", range(prob.shape[1])))
    if 1 in classes:
        idx = classes.index(1)
    elif "1" in classes:
        idx = classes.index("1")
    else:
        idx = prob.shape[1] - 1
    return prob[:, idx].astype(float)


def multiclass_probability_frame(model, X) -> np.ndarray:
    prob = model.predict_proba(X)
    classes = list(getattr(model, "classes_", range(prob.shape[1])))
    aligned = np.zeros((prob.shape[0], len(MULTICLASS_LABELS)), dtype=float)
    if len(classes) == prob.shape[1]:
        for local_idx, cls in enumerate(classes):
            global_idx = int(cls)
            if 0 <= global_idx < len(MULTICLASS_LABELS):
                aligned[:, global_idx] = prob[:, local_idx]
        return aligned
    if prob.shape[1] == len(MULTICLASS_LABELS):
        aligned[:, :] = prob
        return aligned
    raise ValueError(
        f"Cannot align multiclass probabilities with shape {prob.shape} "
        f"to labels {MULTICLASS_LABELS}"
    )

# NANOREMIX_XGB_CLASSES_COMPAT
_nr_load_model_before_classes_fix = load_model

def _nr_repair_xgb_classes(root):
    import numpy as _np
    from xgboost.sklearn import XGBClassifier

    try:
        root_classes = getattr(root, "classes_", None)
    except Exception:
        root_classes = None

    stack = [(root, root_classes)]
    seen = set()

    while stack:
        obj, inherited = stack.pop()
        if obj is None or id(obj) in seen:
            continue
        seen.add(id(obj))

        try:
            own_classes = getattr(obj, "classes_", None)
        except Exception:
            own_classes = None

        classes = own_classes if own_classes is not None else inherited

        if isinstance(obj, XGBClassifier) and own_classes is None:
            if classes is None:
                n_classes = int(getattr(obj, "n_classes_", 2) or 2)
                classes = _np.arange(max(n_classes, 2))
            obj.classes_ = _np.asarray(classes)

        for name in (
            "estimator", "base_estimator", "calibrated_classifiers_",
            "estimators_", "model", "model_", "classifier", "named_steps",
        ):
            try:
                child = getattr(obj, name)
            except Exception:
                continue
            if isinstance(child, dict):
                children = child.values()
            elif isinstance(child, (list, tuple)):
                children = child
            else:
                children = (child,)
            stack.extend((x, classes) for x in children)

    return root

def load_model(path):
    return _nr_repair_xgb_classes(
        _nr_load_model_before_classes_fix(path)
    )
