"""Reference-family conditioning and explicitly separated fusion operations."""

import numpy as np

from .labels import BASE, BINARY, SECOND

EPSILON = 1e-15


def row_validity(p, normalized=True):
    p = np.asarray(p, dtype=float)
    if p.ndim != 2:
        raise ValueError("Expected a two-dimensional probability table")
    valid = np.isfinite(p).all(axis=1) & (p >= 0).all(axis=1) & (p <= 1).all(axis=1)
    if normalized:
        valid &= np.isclose(p.sum(axis=1), 1, atol=1e-6, rtol=0)
    return valid


def logit(p):
    p = np.clip(np.asarray(p, dtype=float), EPSILON, 1 - EPSILON)
    return np.log(p) - np.log1p(-p)


def sigmoid(x):
    x = np.asarray(x, dtype=float)
    return np.exp(-np.logaddexp(0, -x))


def binary_product(probabilities, weights=None, bias=0.0):
    if not probabilities:
        raise ValueError("No probability sources")
    weights = [1.0] * len(probabilities) if weights is None else weights
    if len(weights) != len(probabilities):
        raise ValueError("Fusion weights do not match sources")
    weights = np.asarray(weights, dtype=float)
    if not np.isfinite(weights).all() or (weights < 0).any() or not weights.sum() > 0 or not np.isfinite(bias):
        raise ValueError("Fusion weights must be finite, nonnegative, and nonzero in total")
    z = np.full(np.asarray(probabilities[0]).shape, bias, dtype=float)
    for p, w in zip(probabilities, weights):
        p = np.asarray(p, dtype=float)
        if p.shape != z.shape or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
            raise ValueError("Invalid binary probability values or shape")
        if w:
            z += w * logit(p)
    return sigmoid(z)


def condition(p, source, selected, families):
    p = np.asarray(p, dtype=float)
    q = p[:, [source.index(c) for c in selected]].copy()
    family_mask = np.asarray(families)[:, None] == np.array([BASE[c] for c in selected])[None, :]
    q[~family_mask] = 0.0
    mass = q.sum(axis=1)
    family_mass = np.zeros(len(p))
    for base in "ACGU":
        rows = np.asarray(families) == base
        columns = [i for i, c in enumerate(source) if BASE[c] == base]
        family_mass[rows] = p[np.ix_(rows, columns)].sum(axis=1)
    valid = row_validity(p) & (mass > 0) & (family_mass > 0)
    q[valid] /= mass[valid, None]
    q[~valid] = np.nan
    relative_mass = np.divide(mass, family_mass, out=np.zeros(len(p)), where=family_mass > 0)
    return q, relative_mass, valid


def binary_distribution(pair_probabilities, selected, families):
    q = np.zeros((len(families), len(selected)), dtype=float)
    for base in "ACGU":
        rows = np.asarray(families) == base
        if base in selected:
            q[rows, selected.index(base)] = 1.0
    for mod in BINARY:
        if mod not in selected:
            continue
        if BASE[mod] not in selected:
            raise ValueError("Binary requires the canonical counterpart")
        rows = np.asarray(families) == BASE[mod]
        p = np.asarray(pair_probabilities[mod], dtype=float)
        if p.shape != (len(families),):
            raise ValueError("Binary pair probabilities must retain input row identity")
        q[rows, selected.index(BASE[mod])] = 1 - p[rows]
        q[rows, selected.index(mod)] = p[rows]
    return q, row_validity(q)


def second_pairs(p, selected):
    result = {}
    masses = {}
    for mod in BINARY:
        if mod not in selected:
            continue
        base = BASE[mod]
        mass = p[:, SECOND.index(base)] + p[:, SECOND.index(mod)]
        result[mod] = np.divide(p[:, SECOND.index(mod)], mass, out=np.full(len(p), np.nan), where=mass > 0)
        masses[mod] = mass
    return result, masses


def fuse_multi(third, second, native, selected, families, second_weight=1.0, third_weight=1.0):
    for w in (second_weight, third_weight):
        if not np.isfinite(w) or w <= 0:
            raise ValueError("Multi fusion exponents must be finite and positive")
    r, selected_mass, valid = condition(third, native, selected, families)
    out = r.copy()
    second_mass = np.ones(len(r))
    for base in "ACGU":
        rows = np.flatnonzero((np.asarray(families) == base) & valid)
        common = [c for c in selected if BASE[c] == base and c in SECOND]
        if not len(rows) or not common:
            continue
        columns = [selected.index(c) for c in common]
        shared = r[np.ix_(rows, columns)]
        mass = shared.sum(axis=1)
        active = mass > 0
        rows, shared, mass = rows[active], shared[active], mass[active]
        if not len(rows):
            continue
        q = second[np.ix_(rows, [SECOND.index(c) for c in common])]
        raw_mass = q.sum(axis=1)
        family_columns = [i for i, c in enumerate(SECOND) if BASE[c] == base]
        family_mass = second[np.ix_(rows, family_columns)].sum(axis=1)
        second_mass[rows] = np.divide(raw_mass, family_mass, out=np.zeros(len(rows)), where=family_mass > 0)
        supported = row_validity(second[rows]) & (raw_mass > 0)
        valid[rows[~supported]] = False
        rows, shared, mass, q, raw_mass = rows[supported], shared[supported], mass[supported], q[supported], raw_mass[supported]
        if not len(rows):
            continue
        with np.errstate(divide="ignore"):
            scores = third_weight * np.log(shared / mass[:, None]) + second_weight * np.log(q / raw_mass[:, None])
        maximum = scores.max(axis=1)
        intersects = np.isfinite(maximum)
        valid[rows[~intersects]] = False
        rows, scores, mass, maximum = rows[intersects], scores[intersects], mass[intersects], maximum[intersects]
        combined = np.exp(scores - maximum[:, None])
        combined /= combined.sum(axis=1, keepdims=True)
        out[np.ix_(rows, columns)] = mass[:, None] * combined
    valid &= row_validity(out)
    out[~valid] = np.nan
    return out, selected_mass, second_mass, valid
