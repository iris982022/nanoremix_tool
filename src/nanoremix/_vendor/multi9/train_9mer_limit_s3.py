"""Complete standalone CPU training/inference; load trusted pickles only.

Import this source before direct pickle.load, or use its --predict CLI.
The stable alias below resolves serialized classes without any sibling module.
"""
import argparse
import hashlib
import importlib.metadata
import itertools
import json
import pickle
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import softmax
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import (classification_report, confusion_matrix, log_loss,
                             matthews_corrcoef, roc_auc_score, average_precision_score)
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBClassifier

LABELS5 = 'A C G U m6A m1A Inosine Am m5C hm5C ac4C Cm FormylC pseudoU m1psi m5U Um m7G o8G'.split()
LABELS9 = 'A C G U m6A m1A Am m5C hm5C ac4C o8G Gm m5U pseudoU Um m1psi'.split()
FAMILY = {x: b for b, words in {'A': 'A m6A m1A Inosine Am', 'C': 'C m5C hm5C ac4C Cm FormylC',
          'G': 'G m7G o8G Gm', 'U': 'U pseudoU m1psi m5U Um'}.items() for x in words.split()}
ALIASES = {'T': 'U', **{f'normal{x}': x for x in 'ACGU'}, 'psi': 'pseudoU', 'Psi': 'pseudoU',
           'PSI': 'pseudoU', '\u03a8': 'pseudoU', **{x: 'o8G' for x in ['8-oxoG', '8oxoG', 'oxoG', 'OxoG', 'OG']}}
SPLITS = ['train', 'calibration', 'fusion_validation', 'test']
K, STRATEGY = 9, 's3'
MODULE_ALIAS = 'nr_six_v3_9mer_limit_s3'
STUDY = '9mer_limit_s3'
sys.modules[MODULE_ALIAS] = sys.modules[__name__]
ACTIVE_CLASSES = [c for c in LABELS5 if c != 'Cm'] if K == 5 else LABELS9[:]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=True, allow_nan=False) + '\n')

def token(value):
    value = str(value).strip()
    value = ALIASES.get(value, value)
    if value not in FAMILY:
        raise ValueError(f'Unknown token: {value!r}')
    return value

def schemas(k, strategy):
    neighbors = [i for i in range(1, 2*k) if i != k]
    result = []
    for j in (range(k) if strategy == 's1' else [None]):
        pos = neighbors[j:j+k-1] if j is not None else neighbors
        positions = [f'position{i}' for i in pos] + [f'position{k}_ref']
        numeric = [f'{stat}_output{i}' for i in ([j+1] if j is not None else range(1, k+1))
                   for stat in ['mean', 'median', 'std', 'length']]
        result.append({'name': f'm{j+1}' if j is not None else 'm_all', 'positions': positions,
                       'numeric': numeric, 'features': [f'{p}_{b}' for p in positions for b in 'ACGU'] + numeric})
    return result

def prepare(df, k, classes, label_column=None, training=True):
    df = df.copy()
    if '9mer_final' not in df:
        raise ValueError('Required prepared 9mer_final column absent (9/17 underscore-separated tokens)')
    parts = [str(s).split('_') for s in df['9mer_final']]
    if any(len(p) != 2*k-1 for p in parts):
        raise ValueError(f'9mer_final must contain exactly {2*k-1} tokens')
    annotated = np.array([[token(t) for t in p] for p in parts])
    bases = np.array([[FAMILY[t] for t in p] for p in annotated])
    # Reconstruct; never trust duplicated/mangled position columns from CSV.
    for i in range(2*k-1):
        df[f'position{i+1}'] = bases[:, i]
    df[f'position{k}_ref'] = bases[:, k-1]
    df['canonical_sequence'] = [''.join(p) for p in bases]
    if 'sequence_input' in df:
        for row, (s, b) in enumerate(zip(df.sequence_input, bases[:, k-1])):
            seq = str(s)
            if len(seq) != k or any(t not in 'ACGUT' for t in seq):
                raise ValueError(f'Invalid canonical sequence_input at row {row}')
            if FAMILY[token(seq[k//2])] != b:
                raise ValueError(f'Canonical center conflicts with prepared center at row {row}')
    for info in schemas(k, 's3'):
        for c in info['numeric']:
            if c not in df:
                raise ValueError(f'Missing mandatory numeric column: {c}')
            df[c] = pd.to_numeric(df[c], errors='raise')
            if not np.isfinite(df[c].to_numpy(dtype=float)).all():
                raise ValueError(f'Nonfinite mandatory numeric column: {c}')
    if training:
        if label_column and label_column not in df:
            raise ValueError(f'Missing explicit label column {label_column}')
        labels = np.array([token(t) for t in df[label_column]]) if label_column else annotated[:, k-1]
        if any(FAMILY[t] != b for t, b in zip(labels, bases[:, k-1])):
            raise ValueError('True-label family conflicts with canonical input center')
        if label_column and any(a not in 'ACGU' and a != y for a, y in zip(annotated[:, k-1], labels)):
            raise ValueError('Annotated modification center conflicts with explicit true label')
        unknown = set(labels) - set(classes)
        if unknown:
            raise ValueError(f'Labels outside active namespace: {sorted(unknown)}')
        df['true_label'] = labels
    if 'record_key' not in df:
        if {'src_file', 'src_row'} <= set(df):
            df['record_key'] = [json.dumps([str(a), str(b)]) for a, b in zip(df.src_file, df.src_row)]
        else:
            df['record_key'] = [f'row:{i}' for i in range(len(df))]
    if df.record_key.isna().any() or df.record_key.duplicated().any():
        raise ValueError('record_key must be nonmissing and unique')
    return df

def features(df, info):
    x = np.column_stack([(df[p].to_numpy() == b).astype(float) for p in info['positions'] for b in 'ACGU']
                        + [df[c].to_numpy(dtype=float) for c in info['numeric']])
    return pd.DataFrame(x, columns=info['features'])

def groups_for(df, columns, leakage):
    cols = columns.split(',')
    extra = leakage.split(',') if leakage else []
    for c in cols + extra:
        if c not in df or df[c].isna().any() or (df[c].astype(str).str.strip() == '').any():
            raise ValueError(f'Grouping column missing/empty: {c}')
    keys = [json.dumps(list(row)) for row in df[cols].astype(str).itertuples(index=False, name=None)]
    codes, unique = pd.factorize(np.asarray(keys), sort=True)
    parent = np.arange(len(unique))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for c in extra:
        first = {}
        for i, value in enumerate(df[c].astype(str)):
            r = root(codes[i])
            if value in first:
                a = root(first[value])
                parent[max(a, r)] = min(a, r)
            else:
                first[value] = r
    return np.array([root(c) for c in codes])

def support(y, groups, classes):
    return {c: int(len(np.unique(groups[y == i]))) for i, c in enumerate(classes)}

def allocate(y, groups, classes, folds, seed, budget):
    """Reserve mixed-label group coverage before filling target proportions; validate SGKF."""
    unique, gi = np.unique(groups, return_inverse=True)
    k = len(classes)
    counts = np.zeros((len(unique), k), dtype=int)
    np.add.at(counts, (gi, y), 1)
    present = counts > 0
    sup = present.sum(0)
    if np.any(sup < folds + 3):
        raise ValueError(f'Need at least folds+3={folds+3} groups per class; support={dict(zip(classes, sup.tolist()))}')
    targets = np.array([.55, .15, .15, .15])[:, None] * counts.sum(0)
    for attempt in range(budget):
        rng = np.random.default_rng(seed + attempt)
        assignment = np.full(len(unique), -1)
        coverage = np.zeros((4, k), int)
        totals = np.zeros((4, k), int)
        need = np.tile(np.array([folds, 1, 1, 1])[:, None], (1, k))
        success = True
        while np.any(coverage < need):
            deficits = np.argwhere(coverage < need)
            available = np.array([np.sum((assignment < 0) & present[:, c]) for s, c in deficits])
            if not available.min():
                success = False
                break
            options = np.flatnonzero(available == available.min())
            s, c = deficits[rng.choice(options)]
            choices = np.flatnonzero((assignment < 0) & present[:, c])
            gain = (present[choices] * (coverage[s] < need[s]) / sup).sum(1)
            g = choices[np.argmax(gain + rng.random(len(choices))*1e-6)]
            assignment[g] = s
            coverage[s] += present[g]
            totals[s] += counts[g]
        if not success:
            continue
        for g in rng.permutation(np.flatnonzero(assignment < 0)):
            delta = (((totals + counts[g] - targets)**2 - (totals-targets)**2) / np.maximum(targets, 1)**2).sum(1)
            s = int(np.argmin(delta))
            assignment[g] = s
            totals[s] += counts[g]
        split = assignment[gi]
        tr = np.flatnonzero(split == 0)
        for inner_attempt in range(min(20, budget)):
            cv = list(StratifiedGroupKFold(folds, shuffle=True, random_state=seed+attempt*20+inner_attempt).split(tr, y[tr], groups[tr]))
            if all(len(np.unique(y[tr[a]])) == k and len(np.unique(y[tr[b]])) == k
                   and not set(groups[tr[a]]) & set(groups[tr[b]]) for a, b in cv):
                return split, [(tr[a], tr[b]) for a, b in cv], {'outer_attempt': attempt, 'inner_attempt': inner_attempt}
    raise ValueError(f'Constrained allocation exhausted {budget} outer x <=20 inner attempts; not proof of infeasibility. Change explicit design/budget, not grouping silently.')

def audit(df, y, groups, classes, split=None):
    report = {'rows': len(df), 'groups': len(np.unique(groups)), 'class_group_support': support(y, groups, classes),
              'independence': 'Only specified grouping keys enforced; overlap audit is not proof of biological independence.'}
    if split is not None:
        report['partitions'] = {name: {'rows': int(sum(split == s)), 'class_group_support': support(y[split == s], groups[split == s], classes)} for s, name in enumerate(SPLITS)}
        report['overlap'] = {}
        for col in ['read_id', 'fast5_prefix', 'canonical_sequence']:
            if col not in df:
                continue
            memberships = pd.DataFrame({'value': df[col].astype(str), 'split': split}).drop_duplicates()
            report['overlap'][col] = {f'{SPLITS[a]}:{SPLITS[b]}': len(set(memberships.loc[memberships.split == a, 'value']) & set(memberships.loc[memberships.split == b, 'value'])) for a, b in itertools.combinations(range(4), 2)}
    return report

def train_weights(y, mode='sqrt', beta=.999, supplied=None):
    counts = np.bincount(y)
    n = counts[y].astype(float)
    w = {'none': lambda: np.ones(len(y)), 'sqrt': lambda: 1/np.sqrt(n), 'inverse': lambda: 1/n,
         'effective': lambda: (1-beta)/(-np.expm1(n*np.log(beta)))}[mode]()
    if supplied is not None:
        if mode != 'none':
            raise ValueError('Input weights require --class-weight none to avoid double class weighting')
        w = np.asarray(supplied, float)
    if not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError('Weights must be finite and positive')
    return w / w.mean()

class NamedXGB(ClassifierMixin, BaseEstimator):
    def __init__(self, params=None, n_classes=19):
        self.params = params
        self.n_classes = n_classes

    def fit(self, X, y, sample_weight=None):
        self.classes_ = np.arange(self.n_classes)
        if not np.array_equal(np.unique(y), self.classes_):
            raise ValueError('Missing training classes; no zero-padding or remapping allowed')
        self.model_ = XGBClassifier(**self.params)
        self.model_.fit(X, y, sample_weight=sample_weight)
        self.n_features_in_ = X.shape[1]
        return self

    def predict_proba(self, X):
        if not np.array_equal(self.model_.classes_, self.classes_):
            raise ValueError('Raw model class axis mismatch')
        return validated_probabilities(self.model_.predict_proba(X), self.n_classes)

    def predict(self, X):
        return np.argmax(self.predict_proba(X), axis=1)

def validated_probabilities(p, n_classes):
    """Repair only tiny row-sum drift; never softmax probabilities or reorder columns."""
    p = np.asarray(p, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != n_classes or not len(p):
        raise ValueError('Invalid probability matrix shape')
    if not np.isfinite(p).all() or np.any(p < 0) or np.any(p > 1):
        raise ValueError('Probabilities must be finite and in [0, 1]')
    sums = p.sum(axis=1, keepdims=True)
    if np.any(sums <= 0) or not np.allclose(sums, 1., rtol=0., atol=1e-6):
        raise ValueError('Probability row sums deviate from one by more than 1e-6')
    return p / sums

def probability_log_loss(y, p, classes):
    p = validated_probabilities(p, len(classes))
    # sklearn 1.3 checks sums AFTER clipping: gated zeros add K * epsilon.
    # Normalize a metric-only clipped copy; keep exact zeros in inference.
    eps = np.finfo(np.float64).eps
    p = np.clip(p, eps, 1 - eps)
    p /= p.sum(axis=1, keepdims=True)
    return log_loss(y, p, labels=np.arange(len(classes)))

def group_nll(y, p, groups):
    losses = -np.log(np.clip(p[np.arange(len(y)), y], 1e-15, 1))
    return float(pd.DataFrame({'g': groups, 'loss': losses}).groupby('g').loss.mean().mean())

def aggregate(probs, bases, classes, config):
    m, n, k = probs.shape
    out = np.zeros((n, k))
    for b in 'ACGU':
        rows = np.flatnonzero(bases == b)
        ids = np.flatnonzero([FAMILY[c] == b for c in classes]) if config['gate'] else np.arange(k)
        if not len(rows):
            continue
        if not len(ids):
            raise ValueError(f'No active labels for canonical input family {b}')
        p = np.clip(probs[:, rows][:, :, ids], 1e-15, 1)
        p /= p.sum(2, keepdims=True)
        spec = config['families'][b]
        w = np.asarray(spec['weights'])
        if spec['method'] == 'arithmetic':
            z = np.log(np.einsum('m,mnk->nk', w, p))
        else:
            z = np.einsum('m,mnk->nk', w, np.log(p))
        out[np.ix_(rows, ids)] = softmax(z/spec['temperature'], axis=1)
    return out

def fixed_config(m, method='geometric', temperature=1., gate=True):
    return {'gate': gate, 'families': {b: {'method': method, 'weights': [1/m]*m, 'temperature': temperature} for b in 'ACGU'}}

def tune_fusion(probs, bases, y, groups, classes, strategy, penalty):
    m = len(probs)
    candidates = {'original_product': fixed_config(m, temperature=1/m), 'arithmetic': fixed_config(m, 'arithmetic'),
                  'geometric': fixed_config(m)}
    optimized = fixed_config(m)
    diagnostics = {}
    for b in 'ACGU':
        rows = np.flatnonzero(bases == b)
        def objective(x):
            spec = optimized['families'][b]
            spec['weights'] = (x[:-1] if strategy == 's1' else np.ones(1)).tolist()
            spec['temperature'] = float(x[-1])
            p = aggregate(probs[:, rows], bases[rows], classes, optimized)
            reg = penalty * (np.sum((np.asarray(spec['weights'])-1/m)**2) + np.log(x[-1])**2)
            return group_nll(y[rows], p, groups[rows]) + reg
        x0 = np.r_[np.ones(m)/m, 1.] if strategy == 's1' else np.ones(1)
        result = minimize(objective, x0, method='SLSQP', bounds=([(0, 1)]*m if strategy == 's1' else []) + [(.25, 4)],
                          constraints=[{'type': 'eq', 'fun': lambda x: x[:-1].sum()-1}] if strategy == 's1' else [],
                          options={'maxiter': 150, 'ftol': 1e-9})
        if not result.success:
            raise ValueError(f'Fusion optimizer failed for {b}: {result.message}')
        objective(result.x)
        diagnostics[b] = {'objective_with_penalty': float(result.fun), 'iterations': int(result.nit)}
    candidates['regularized'] = optimized
    scores = {name: group_nll(y, aggregate(probs, bases, classes, cfg), groups) for name, cfg in candidates.items()}
    winner = min(scores, key=scores.get)
    return candidates[winner], {'selection': winner, 'group_nll': scores, 'candidates': candidates, 'optimizer': diagnostics,
                                'source_partition': 'fusion_validation', 'penalty': penalty}

class ModelBundle:
    def __init__(self, manifest, components):
        self.manifest = manifest
        self.components = components
        self.classes_ = np.array(manifest['classes'])

    def predict_proba(self, df):
        d = prepare(df, self.manifest['k'], self.manifest['classes'], training=False)
        probs = np.array([m.predict_proba(features(d, s)) for m, s in zip(self.components, self.manifest['schemas'])])
        return aggregate(probs, d[f'position{self.manifest["k"]}_ref'].to_numpy(), self.manifest['classes'], self.manifest['fusion'])

    def predict(self, df):
        return self.classes_[self.predict_proba(df).argmax(1)]

def load_model(directory):
    directory = Path(directory)
    manifest = json.loads((directory/'manifest.json').read_text())
    if (manifest['format'] != 'nanoremix-thirdgen-standalone-v1'
            or manifest['module_alias'] != MODULE_ALIAS
            or manifest['k'] != K or manifest['strategy'] != STRATEGY
            or manifest['classes'] != ACTIVE_CLASSES):
        raise ValueError('Model must match this standalone script and fixed namespace')
    models = []
    for item in manifest['components']:
        path = directory/item['pickle']
        if sha(path) != item['sha256']:
            raise ValueError(f'Model hash mismatch: {path}')
        with path.open('rb') as f:
            model = pickle.load(f)
        if not np.array_equal(model.classes_, np.arange(len(manifest['classes']))):
            raise ValueError('Component class axis mismatch')
        models.append(model)
    return ModelBundle(manifest, models)

def evaluate(y, p, classes, bases, groups, out, name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pred = p.argmax(1)
    result = {'classification': classification_report(y, pred, labels=np.arange(len(classes)), target_names=classes,
               output_dict=True, zero_division=0), 'log_loss': float(probability_log_loss(y, p, classes)),
              'group_nll': group_nll(y, p, groups), 'MCC': float(matthews_corrcoef(y, pred)),
              'Brier_multiclass_sum': float(np.mean(np.sum((p-np.eye(len(classes))[y])**2, axis=1)))}
    ranking = {}
    for family in ['overall', *'ACGU']:
        rows = np.arange(len(y)) if family == 'overall' else np.flatnonzero(bases == family)
        ids = list(range(len(classes))) if family == 'overall' else [i for i, c in enumerate(classes) if FAMILY[c] == family]
        values = {}
        q = p[np.ix_(rows, ids)]
        q = q/np.maximum(q.sum(1, keepdims=True), 1e-15)
        for j, c in enumerate(ids):
            truth = y[rows] == c
            values[classes[c]] = {'AUC': float(roc_auc_score(truth, q[:, j])), 'AP': float(average_precision_score(truth, q[:, j]))} if len(np.unique(truth)) == 2 else {'AUC': None, 'AP': None}
        valid = [v for v in values.values() if v['AUC'] is not None]
        ranking[family] = {'per_class': values, 'macro_AUC': float(np.mean([v['AUC'] for v in valid])) if valid else None,
                           'macro_AP': float(np.mean([v['AP'] for v in valid])) if valid else None}
    result['ranking'] = ranking
    dump_json(out/f'metrics_{name}.json', result)
    cm = confusion_matrix(y, pred, labels=np.arange(len(classes)))
    for suffix, matrix in [('count', cm), ('rowpercent', cm/np.maximum(cm.sum(1, keepdims=True), 1)*100)]:
        pd.DataFrame(matrix, index=classes, columns=classes).to_csv(out/f'confusion_{name}_{suffix}.csv')
        fig, ax = plt.subplots(figsize=(11, 9))
        im = ax.imshow(matrix, cmap='Blues')
        ax.set(xticks=range(len(classes)), yticks=range(len(classes)), xticklabels=classes, yticklabels=classes,
               xlabel='Predicted', ylabel='True', title=f'{name} {suffix}')
        plt.setp(ax.get_xticklabels(), rotation=90)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(out/f'confusion_{name}_{suffix}.png', dpi=120)
        plt.close(fig)
    return result

def reuse_partitions(path, d, y, groups, manifest, folds):
    path = Path(path)
    baseline = json.loads(path.with_name('run_config.json').read_text())
    for key in ['input_sha256', 'k', 'strategy', 'classes', 'group_semantics']:
        if baseline[key] != manifest[key]:
            raise ValueError(f'Baseline mismatch: {key}')
    for key in ['group_columns', 'leakage_columns', 'label_column', 'folds', 'seed',
                'class_weight', 'beta', 'input_weight_column', 'calibration_weight_column']:
        if baseline['args'].get(key) != manifest['args'].get(key):
            raise ValueError(f'Baseline argument mismatch: {key}')
    p = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = ['record_key', 'true_label', 'canonical_sequence', 'effective_group',
                'partition', 'inner_validation_fold'] + manifest['args']['group_columns'].split(',')
    if not set(required) <= set(p) or p.record_key.duplicated().any():
        raise ValueError('Missing partition columns or duplicate record keys')
    keys = d.record_key.astype(str)
    if keys.duplicated().any() or set(keys) != set(p.record_key):
        raise ValueError('Partition record keys do not match input exactly')
    p = p.set_index('record_key').loc[keys].reset_index()
    for col in ['true_label', 'canonical_sequence'] + manifest['args']['group_columns'].split(','):
        if not np.array_equal(p[col].to_numpy(), d[col].astype(str).to_numpy()):
            raise ValueError(f'Partition input mismatch: {col}')
    if not np.array_equal(p.effective_group.to_numpy(), groups.astype(str)):
        raise ValueError('Effective grouping mismatch')
    if not p.partition.isin(SPLITS).all():
        raise ValueError('Unknown partition membership')
    split = np.array([SPLITS.index(s) for s in p.partition])
    fold = pd.to_numeric(p.inner_validation_fold, errors='raise').to_numpy()
    if not np.isin(fold[split == 0], np.arange(folds)).all() or np.any(fold[split != 0] != -1):
        raise ValueError('Invalid inner fold membership')
    for membership in [split, np.where(split == 0, fold, -1)]:
        if pd.DataFrame({'g': groups, 's': membership}).groupby('g').s.nunique().max() != 1:
            raise ValueError('Group crosses partition or inner fold')
    classes = np.arange(len(manifest['classes']))
    for s in range(4):
        if not np.array_equal(np.unique(y[split == s]), classes):
            raise ValueError('Missing class coverage in fixed partition')
    cv = [(np.flatnonzero((split == 0) & (fold != f)),
           np.flatnonzero((split == 0) & (fold == f))) for f in range(folds)]
    for a, b in cv:
        if any(not np.array_equal(np.unique(y[ids]), classes) for ids in [a, b]):
            raise ValueError('Missing class coverage in fixed inner fold')
    return split, cv, {'reused': True, 'partitions_file': str(path.resolve()),
                       'partitions_sha256': sha(path),
                       'baseline_config_sha256': sha(path.with_name('run_config.json'))}

def main(k, strategy):
    parser = argparse.ArgumentParser(description=f'{k}mer {strategy}: independent CPU retraining')
    parser.add_argument('--input')
    parser.add_argument('--output-dir', '--output_dir', dest='output_dir')
    parser.add_argument('--group-columns', default='fast5_prefix' if k == 5 else 'read_id',
                        help='Default: fast5_prefix (5mer), read_id (9mer). Tuple keys; canonical overlap audited, not excluded.')
    parser.add_argument('--leakage-columns', default='', help='Explicit extra columns joined by transitive union only')
    parser.add_argument('--label-column', help='Explicit named true label; never infer numeric label_id')
    parser.add_argument('--predict', action='store_true', help='Inference only; requires --model-dir and --output')
    parser.add_argument('--model-dir', help='Trusted model directory from this exact standalone variant')
    parser.add_argument('--output', help='New prediction CSV (never overwritten)')
    parser.add_argument('--n-trials', type=int, default=40)
    parser.add_argument('--max-trees', type=int, default=300)
    parser.add_argument('--depth-min', type=int, default=3)
    parser.add_argument('--depth-max', type=int, default=10)
    parser.add_argument('--tree-grid', default='300,600', help='Predeclared tree counts; overrides --max-trees (use a singleton for fixed trees)')
    parser.add_argument('--partitions-file', help='Existing partitions.csv with adjacent baseline run_config.json')
    parser.add_argument('--tuning-only', action='store_true', help='Fixed-partition Train CV only; no final model or test scoring')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--folds', type=int, default=3)
    parser.add_argument('--seed', type=int, default=25)
    parser.add_argument('--split-budget', type=int, default=100)
    parser.add_argument('--class-weight', choices=['none', 'sqrt', 'inverse', 'effective'], default='sqrt')
    parser.add_argument('--beta', type=float)
    parser.add_argument('--input-weight-column', help='Opt-in train-only weights; requires class-weight none')
    parser.add_argument('--calibration-weight-column', help='Explicit deployment weights, never automatic class balance')
    parser.add_argument('--fusion-penalty', type=float, default=.02)
    parser.add_argument('--save-raw-probs', action='store_true')
    parser.add_argument('--check-data', action='store_true')
    parser.add_argument('--check-dependencies', action='store_true')
    args = parser.parse_args()
    env = {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'xgboost', 'optuna', 'matplotlib']}
    env.update(python=sys.version, executable=sys.executable, platform=platform.platform())
    if args.check_dependencies:
        print(json.dumps(env, indent=2)); return
    if args.predict or args.model_dir:
        if not args.input or not args.model_dir or not args.output:
            parser.error('Inference requires --input, --model-dir and --output')
        destination = Path(args.output)
        if destination.exists():
            parser.error('Refusing to overwrite prediction output')
        model = load_model(args.model_dir)
        d = prepare(pd.read_csv(args.input, low_memory=False), k, ACTIVE_CLASSES, training=False)
        probs = model.predict_proba(d)
        result = pd.DataFrame(probs, columns=[f'prob_{c}' for c in model.classes_])
        result.insert(0, 'record_key', d.record_key.to_numpy())
        result['prediction'] = model.classes_[probs.argmax(1)]
        with destination.open('x') as f:
            result.to_csv(f, index=False)
        print(f'Predictions: {destination}')
        return
    if not args.input or not args.output_dir or not args.group_columns:
        parser.error('--input, --output-dir and nonempty --group-columns are required')
    if min(args.n_trials, args.max_trees, args.threads, args.split_budget) < 1 or args.folds < 2 or args.fusion_penalty < 0:
        parser.error('Invalid training budget/folds/penalty')
    if not 1 <= args.depth_min <= args.depth_max:
        parser.error('Require 1 <= depth-min <= depth-max')
    tree_grid = None
    if args.tree_grid is not None:
        try:
            tree_grid = [int(s) for s in args.tree_grid.split(',')]
            if min(tree_grid) < 1 or len(set(tree_grid)) != len(tree_grid):
                raise ValueError()
        except ValueError:
            parser.error('--tree-grid requires unique positive integers')
    if args.tuning_only and (not args.partitions_file or args.save_raw_probs):
        parser.error('--tuning-only requires --partitions-file and forbids --save-raw-probs')
    if args.beta is not None and args.class_weight != 'effective':
        warnings.warn('--beta ignored unless --class-weight effective')
    beta = .999 if args.beta is None else args.beta
    if args.class_weight == 'effective' and not 0 < beta < 1:
        parser.error('effective beta must lie in (0,1)')
    if args.input_weight_column and args.class_weight != 'none':
        parser.error('--input-weight-column requires --class-weight none')
    out = Path(args.output_dir).resolve()
    protected = [Path(args.input).resolve().parent]
    if args.partitions_file:
        protected.append(Path(args.partitions_file).resolve().parent)
    if out.exists() or any(out == p or p in out.parents or out in p.parents for p in protected):
        parser.error('Output must not exist and must be outside input/baseline directories')
    out.mkdir(parents=True, exist_ok=False)
    classes0 = LABELS5 if k == 5 else LABELS9
    excluded = ['Cm'] if k == 5 else []
    classes = ACTIVE_CLASSES[:]
    print(f'Fixed configuration: {k}mer {strategy}, {len(classes)} classes, excluded={excluded}, groups={args.group_columns}', flush=True)
    manifest = {'format': 'nanoremix-thirdgen-standalone-v1', 'module_alias': MODULE_ALIAS,
                'k': k, 'strategy': strategy, 'classes': classes, 'study': STUDY,
                'evaluation_context': 'historical-test development comparison; not new external validation',
                'namespace': f'{k}mer_fixed{len(classes)}',
                'excluded_labels': excluded, 'args': vars(args), 'environment': env,
                'input_sha256': sha(args.input), 'schemas': schemas(k, strategy),
                'scripts_sha256': {Path(__file__).name: sha(__file__)},
                'calibration': 'sigmoid; calibration partition only; unweighted unless explicit deployment column',
                'weighting': 'fold training labels only; row mean 1; no oversampling',
                'group_semantics': 'tuple of group-columns; additional union only explicit leakage-columns'}
    dump_json(out/'run_config.json', manifest)
    try:
        d = prepare(pd.read_csv(args.input, low_memory=False), k, classes0, args.label_column)
        excluded_rows = int(d.true_label.isin(excluded).sum())
        d = d.loc[~d.true_label.isin(excluded)].reset_index(drop=True)
        manifest['excluded_rows'] = excluded_rows
        y = np.array([classes.index(c) for c in d.true_label])
        groups = groups_for(d, args.group_columns, args.leakage_columns)
        dump_json(out/'group_audit.json', audit(d, y, groups, classes))
        for col in [args.input_weight_column, args.calibration_weight_column]:
            if col:
                train_weights(y, 'none', supplied=pd.to_numeric(d[col], errors='raise'))
        if args.partitions_file:
            split, cv, search = reuse_partitions(args.partitions_file, d, y, groups, manifest, args.folds)
        else:
            split, cv, search = allocate(y, groups, classes, args.folds, args.seed, args.split_budget)
        manifest['split_search'] = search
        dump_json(out/'run_config.json', manifest)
        dump_json(out/'group_audit.json', audit(d, y, groups, classes, split))
        partition = d[['record_key', 'true_label', 'canonical_sequence']].copy()
        for c in dict.fromkeys(args.group_columns.split(',') + ['read_id', 'fast5_prefix']):
            if c in d:
                partition[c] = d[c]
        partition['effective_group'] = groups
        partition['partition'] = np.array(SPLITS)[split]
        partition['inner_validation_fold'] = -1
        for f, (a, b) in enumerate(cv):
            partition.loc[b, 'inner_validation_fold'] = f
        partition.to_csv(out/'partitions.csv', index=False)
        dump_json(out/'feasibility.json', {'feasible': True, 'search': search, 'classes': classes, 'excluded_rows': excluded_rows})
    except (ValueError, KeyError) as e:
        dump_json(out/'feasibility.json', {'feasible': False, 'error': str(e), 'classes': classes})
        raise
    if args.check_data:
        print(f'Feasible: {out}'); return
    import optuna
    tr, cal, val, test = [np.flatnonzero(split == s) for s in range(4)]
    models, validation_probs, components = [], [], []
    def weights(ids):
        return train_weights(y[ids], args.class_weight, beta, d.iloc[ids][args.input_weight_column] if args.input_weight_column else None)
    for info in manifest['schemas']:
        x = features(d, info)
        fixed = {'objective': 'multi:softprob', 'num_class': len(classes), 'eval_metric': 'mlogloss',
                 'tree_method': 'hist', 'device': 'cpu', 'n_jobs': args.threads, 'random_state': args.seed,
                 'n_estimators': args.max_trees, 'verbosity': 0}
        def objective(trial):
            params = dict(fixed, max_depth=trial.suggest_int('max_depth', args.depth_min, args.depth_max),
                          learning_rate=trial.suggest_float('learning_rate', .02, .2, log=True),
                          min_child_weight=trial.suggest_float('min_child_weight', 1, 8),
                          subsample=trial.suggest_float('subsample', .7, 1),
                          colsample_bytree=trial.suggest_float('colsample_bytree', .7, 1),
                          reg_lambda=trial.suggest_float('reg_lambda', .1, 10, log=True), reg_alpha=0., gamma=0.)
            if tree_grid is not None:
                params['n_estimators'] = trial.suggest_categorical('n_estimators', tree_grid)
            trial.set_user_attr('full_params', params)
            losses = []
            for a, b in cv:
                model = NamedXGB(params, len(classes)).fit(x.iloc[a], y[a], weights(a))
                losses.append(group_nll(y[b], model.predict_proba(x.iloc[b]), groups[b]))
            trial.set_user_attr('fold_losses', losses)
            return float(np.mean(losses))
        study = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler(seed=args.seed))
        study.optimize(objective, n_trials=args.n_trials, n_jobs=1)
        params = study.best_trial.user_attrs['full_params']
        study.trials_dataframe().to_csv(out/f'{info["name"]}_trials.csv', index=False)
        dump_json(out/f'{info["name"]}_best.json', {
            'schema': info, 'params': params, 'trial': study.best_trial.number,
            'fold_losses': study.best_trial.user_attrs['fold_losses'],
            'objective': 'mean_inner_fold_group_nll', 'value': study.best_value,
            'selection_partition': 'train', 'warm_start': False})
        if args.tuning_only:
            continue
        raw = NamedXGB(params, len(classes)).fit(x.iloc[tr], y[tr], weights(tr))
        raw.model_.save_model(out/f'{info["name"]}_raw_UNCALIBRATED.ubj')
        calibrated = CalibratedClassifierCV(raw, method='sigmoid', cv='prefit')
        calibrated.fit(x.iloc[cal], y[cal], sample_weight=d.iloc[cal][args.calibration_weight_column].to_numpy(float) if args.calibration_weight_column else None)
        filename = f'{info["name"]}_calibrated.pkl'
        with (out/filename).open('wb') as f:
            pickle.dump(calibrated, f, protocol=5)
        models.append(calibrated)
        components.append({'pickle': filename, 'sha256': sha(out/filename), 'params': params,
                           'raw_ubj': f'{info["name"]}_raw_UNCALIBRATED.ubj', 'classes': classes})
        validation_probs.append(calibrated.predict_proba(x.iloc[val]))
    if args.tuning_only:
        dump_json(out/'TUNING_COMPLETE.json', {
            'tuning_complete': True, 'completed_model': False, 'test_scored': False,
            'calibration_fitted': False, 'fusion_fitted': False,
            'selection_partition': 'train', 'components': [s['name'] for s in manifest['schemas']],
            'run_config_sha256': sha(out/'run_config.json')})
        print(f'Tuning only completed; no final model or test scores: {out}')
        return
    bases = d[f'position{k}_ref'].to_numpy()
    fusion, selection = tune_fusion(np.array(validation_probs), bases[val], y[val], groups[val], classes, strategy, args.fusion_penalty)
    manifest.update(fusion=fusion, fusion_selection=selection, components=components, frozen_before_test=True)
    dump_json(out/'manifest.json', manifest)
    bundle = ModelBundle(manifest, models)
    with (out/'final_bundle.pkl').open('wb') as f:
        pickle.dump(bundle, f, protocol=5)
    # Test probabilities computed once, only after configuration is persisted and frozen.
    probs = np.array([model.predict_proba(features(d.iloc[test], info)) for model, info in zip(models, manifest['schemas'])])
    configs = {'selected': fusion, 'original_ungated': fixed_config(len(models), temperature=1/len(models), gate=False),
               'original_hardgate': fixed_config(len(models), temperature=1/len(models)),
               'arithmetic_hardgate': fixed_config(len(models), 'arithmetic')}
    for name, cfg in configs.items():
        p = aggregate(probs, bases[test], classes, cfg)
        evaluate(y[test], p, classes, bases[test], groups[test], out, name)
        if name == 'selected':
            pred = pd.DataFrame(p, columns=[f'prob_{c}' for c in classes])
            pred.insert(0, 'record_key', d.iloc[test].record_key.to_numpy())
            pred['true_label'] = d.iloc[test].true_label.to_numpy()
            pred['prediction'] = np.array(classes)[p.argmax(1)]
            pred.to_csv(out/'test_predictions.csv', index=False)
    if args.save_raw_probs:
        np.savez_compressed(out/'internal_probabilities.npz', validation=np.array(validation_probs), test=probs,
                            validation_keys=d.iloc[val].record_key.to_numpy(str), test_keys=d.iloc[test].record_key.to_numpy(str))
    dump_json(out/'completed.json', {'test_evaluated_once_after_freeze': True, 'manifest_sha256': sha(out/'manifest.json')})
    print(f'Completed: {out}')

# Both custom classes can occur in pickles; never serialize them as __main__.
NamedXGB.__module__ = MODULE_ALIAS
ModelBundle.__module__ = MODULE_ALIAS

if __name__ == '__main__':
    main(K, STRATEGY)
