from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

from .constants import DEFAULT_MODEL_ROOT, MULTICLASS_LABELS

try:
    from sklearn.base import BaseEstimator, ClassifierMixin
except Exception:
    BaseEstimator = object
    ClassifierMixin = object


class LocalXGB(ClassifierMixin, BaseEstimator):
    """Compatibility shim for multiclass pkls saved from __main__.LocalXGB."""

    _estimator_type = "classifier"

    def __init__(self, base_params=None, present_classes=None, num_global=None):
        self.base_params = dict(base_params or {})
        present = [] if present_classes is None else present_classes
        self.present = np.array(sorted(int(c) for c in np.unique(present)))
        self.g2l = {g: i for i, g in enumerate(self.present)}
        self.l2g = {i: g for i, g in enumerate(self.present)}
        self.num_global = int(num_global or len(self.present) or len(MULTICLASS_LABELS))
        self.base_params.update({"objective": "multi:softprob", "num_class": len(self.present)})
        self.model = None
        self.classes_ = np.arange(self.num_global)

    def fit(self, X, y, sample_weight=None, verbose=False):
        from xgboost import XGBClassifier

        y_local = np.array([self.g2l[int(c)] for c in y])
        self.model = XGBClassifier(**self.base_params)
        self.model.fit(X, y_local, sample_weight=sample_weight, verbose=verbose)
        return self

    def predict_proba(self, X):
        if self.model is None:
            raise ValueError("LocalXGB.model is not loaded")
        _patch_xgb_pickle_compat(self)
        p_loc = self.model.predict_proba(X)
        P = np.zeros((p_loc.shape[0], self.num_global), dtype=float)
        for local_idx, global_idx in self.l2g.items():
            P[:, global_idx] = p_loc[:, local_idx]
        return P

    def predict(self, X):
        return np.argmax(self.predict_proba(X), axis=1)


def install_pickle_shims() -> None:
    main_mod = sys.modules.get("__main__")
    if main_mod is not None and not hasattr(main_mod, "LocalXGB"):
        setattr(main_mod, "LocalXGB", LocalXGB)


def load_model(path: str | Path):
    install_pickle_shims()
    try:
        import joblib
    except ImportError as exc:
        raise RuntimeError(
            "joblib is required to load NanoRemix model pkls. Install with `pip install -e .`."
        ) from exc
    return joblib.load(path)


def normalize_model_root(model_root: str | Path | None) -> Path:
    raw = model_root or os.environ.get("NANOREMIX5MER_MODEL_ROOT") or DEFAULT_MODEL_ROOT
    root = Path(raw).expanduser().resolve()
    if (root / "binary_classify").is_dir() and (
        (root / "mutli_classify").is_dir() or (root / "multi_classify").is_dir()
    ):
        return root
    if root.name in {"binary_classify", "mutli_classify", "multi_classify"}:
        return root.parent
    if (root / "5mer" / "binary_classify").is_dir():
        return root / "5mer"
    raise FileNotFoundError(
        f"Cannot locate NanoRemix 5mer model root from {root}. "
        "Pass --model-root to the directory containing binary_classify and mutli_classify."
    )


def _find_strategy_dir(parent: Path, strategy: str) -> Path:
    candidate = parent / strategy
    if candidate.is_dir():
        return candidate
    for child in parent.iterdir():
        if child.is_dir() and child.name.strip() == strategy:
            return child
    raise FileNotFoundError(f"Strategy directory not found: {parent}/{strategy}")


def resolve_model_paths(
    model_root: str | Path | None,
    mode: str,
    strategy: str,
    condition: str = "full",
    pair: str | None = None,
) -> dict[str, Path]:
    root = normalize_model_root(model_root)
    if mode == "binary":
        if pair is None:
            raise ValueError("--pair is required for binary mode")
        base_dir = root / "binary_classify" / pair
        strategy_dir = _find_strategy_dir(base_dir, strategy)
        model_dir = strategy_dir if condition == "full" else strategy_dir / f"only_{condition}"
    elif mode == "multiclass":
        if condition != "full":
            raise ValueError("Multiclass pkls are fixed full-feature models; use --condition full.")
        multi_root = root / "mutli_classify"
        if not multi_root.is_dir():
            multi_root = root / "multi_classify"
        nested = multi_root / "model" / strategy
        expected = ["m_all"] if strategy == "strategy3" else [f"m{i}" for i in range(1, 6)]
        model_dir = (
            nested
            if nested.is_dir()
            and all((nested / f"{name}_calibrated.pkl").is_file() for name in expected)
            else multi_root / strategy
        )
    else:
        raise ValueError(f"Unknown mode: {mode}")

    model_names = ["m_all"] if strategy == "strategy3" else [f"m{i}" for i in range(1, 6)]
    paths = {name: model_dir / f"{name}_calibrated.pkl" for name in model_names}
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
        estimator = getattr(calibrated[0], "estimator", None) or getattr(calibrated[0], "base_estimator", None)
        names = getattr(estimator, "feature_names_in_", None)
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




_XGB_COMPAT_DEFAULTS = {
    "use_label_encoder": False,
    "gpu_id": None,
    "predictor": None,
    "feature_types": None,
    "enable_categorical": False,
    "max_cat_to_onehot": None,
    "max_cat_threshold": None,
    "callbacks": None,
    "early_stopping_rounds": None,
    "eval_metric": None,
    "validate_parameters": None,
    "n_jobs": None,
    "missing": float("nan"),
    "importance_type": None,
    "base_score": None,
    "booster": None,
    "objective": None,
    "tree_method": None,
    "random_state": None,
}

def _ensure_xgb_attrs(obj):
    try:
        mod = obj.__class__.__module__.lower()
        name = obj.__class__.__name__
        is_xgb = "xgboost" in mod or name.startswith("XGB")
    except Exception:
        is_xgb = False

    if not is_xgb:
        return

    try:
        import numpy as _np
        if not hasattr(obj, "classes_"):
            n_classes = getattr(obj, "n_classes_", None)
            try:
                n_classes = int(n_classes)
            except Exception:
                n_classes = 2
            if n_classes <= 0:
                n_classes = 2
            obj.classes_ = _np.arange(n_classes)
        if not hasattr(obj, "n_classes_"):
            try:
                obj.n_classes_ = len(obj.classes_)
            except Exception:
                obj.n_classes_ = 2
    except Exception:
        pass

    names = set(_XGB_COMPAT_DEFAULTS)
    try:
        names.update(obj.__class__._get_param_names())
    except Exception:
        pass

    for key in names:
        if not hasattr(obj, key):
            try:
                setattr(obj, key, _XGB_COMPAT_DEFAULTS.get(key, None))
            except Exception:
                pass

def _patch_xgb_pickle_compat(model):
    """Patch old pickled XGBoost sklearn estimators for newer xgboost versions."""
    try:
        import xgboost.sklearn as xgbsk
    except Exception:
        xgbsk = None

    if xgbsk is not None:
        try:
            if not getattr(xgbsk.XGBModel.get_params, "_nanoremix_xgb_compat2", False):
                def get_params_compat(self, deep=True):
                    _ensure_xgb_attrs(self)
                    params = {}

                    try:
                        names = self.__class__._get_param_names()
                    except Exception:
                        names = list(_XGB_COMPAT_DEFAULTS)

                    for key in names:
                        if not hasattr(self, key):
                            try:
                                setattr(self, key, _XGB_COMPAT_DEFAULTS.get(key, None))
                            except Exception:
                                pass

                        try:
                            value = getattr(self, key)
                        except Exception:
                            value = _XGB_COMPAT_DEFAULTS.get(key, None)

                        params[key] = value

                        if deep and hasattr(value, "get_params"):
                            try:
                                for sub_key, sub_val in value.get_params().items():
                                    params[f"{key}__{sub_key}"] = sub_val
                            except Exception:
                                pass

                    try:
                        kwargs = getattr(self, "kwargs", None)
                        if isinstance(kwargs, dict):
                            params.update(kwargs)
                    except Exception:
                        pass

                    return params

                get_params_compat._nanoremix_xgb_compat2 = True
                xgbsk.XGBModel.get_params = get_params_compat
        except Exception:
            pass

        try:
            if not getattr(xgbsk.XGBModel._can_use_inplace_predict, "_nanoremix_xgb_compat2", False):
                def no_inplace_predict(self):
                    _ensure_xgb_attrs(self)
                    return False
                no_inplace_predict._nanoremix_xgb_compat2 = True
                xgbsk.XGBModel._can_use_inplace_predict = no_inplace_predict
        except Exception:
            pass

    seen = set()

    def walk(obj):
        oid = id(obj)
        if oid in seen:
            return
        seen.add(oid)

        _ensure_xgb_attrs(obj)

        for attr in (
            "estimator", "estimator_", "base_estimator", "base_estimator_",
            "_estimator", "classifier", "classifier_", "calibrated_classifiers_",
            "calibrators", "calibrators_", "calibrated_classifiers",
        ):
            try:
                child = getattr(obj, attr)
            except Exception:
                continue

            if isinstance(child, (list, tuple)):
                for x in child:
                    walk(x)
            elif isinstance(child, dict):
                for x in child.values():
                    walk(x)
            else:
                walk(child)

    walk(model)


def binary_positive_probability(model, X) -> np.ndarray:
    _patch_xgb_pickle_compat(model)
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
    _patch_xgb_pickle_compat(model)
    prob = model.predict_proba(X)
    classes = list(getattr(model, "classes_", range(prob.shape[1])))
    aligned = np.zeros((prob.shape[0], len(MULTICLASS_LABELS)), dtype=float)
    if prob.shape[1] == len(MULTICLASS_LABELS) and len(classes) != prob.shape[1]:
        aligned[:, :] = prob
        return aligned
    for local_idx, cls in enumerate(classes):
        global_idx = int(cls)
        if 0 <= global_idx < len(MULTICLASS_LABELS):
            aligned[:, global_idx] = prob[:, local_idx]
    return aligned
