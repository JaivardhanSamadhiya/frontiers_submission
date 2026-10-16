"""Retrospective screening diagnostics; no fitting or threshold optimization."""
from math import comb
import numpy as np


def screening_curve(labels, scores, identifiers):
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    ids = np.asarray(identifiers, dtype=str)
    if not len(y) or not (len(y) == len(s) == len(ids)):
        raise ValueError('Aligned nonempty arrays required')
    if not np.isin(y, [0, 1]).all() or not np.isfinite(s).all():
        raise ValueError('Observed binary outcomes and finite scores required')
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate candidate identifiers')
    order = np.lexsort((ids, -s))
    found = np.cumsum(y[order])
    n, positive = len(y), int(y.sum())
    budget = np.arange(1, n + 1)
    hit = (found > 0).astype(float)
    return {
        'order': ids[order], 'hit': hit, 'precision': found / budget,
        'recall': found / positive if positive else np.full(n, np.nan),
        'first_hit': float(np.flatnonzero(hit)[0] + 1) if positive else np.nan,
        'normalized_hit_area': float(hit.mean()),
        'positive': positive, 'candidates': n,
    }


def random_curve(n, positive):
    if n < 1 or not 0 <= positive <= n:
        raise ValueError('Invalid candidate/positive counts')
    hit = np.array([1 - (comb(n-positive, k) / comb(n, k) if k <= n-positive else 0)
                    for k in range(1, n+1)])
    return {'hit': hit, 'precision': np.full(n, positive/n),
            'recall': np.arange(1, n+1)/n if positive else np.full(n, np.nan),
            'first_hit': (n+1)/(positive+1) if positive else np.nan,
            'normalized_hit_area': float(hit.mean())}


def cluster_mean_ci(values, groups, repetitions=2000, seed=20261005):
    """Percentile cluster bootstrap, row-weighted, missing outcomes excluded.

    For paired differences pass the per-row difference, not separate bootstrap
    intervals. Clusters are resampled within the supplied analysis population.
    """
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    valid = np.isfinite(values)
    if not valid.any():
        return {'mean': None, 'lo': None, 'hi': None, 'n': 0, 'groups': 0,
                'valid_replicates': 0}
    values, groups = values[valid], groups[valid]
    unique, inverse = np.unique(groups, return_inverse=True)
    result = {'mean': float(values.mean()), 'n': len(values), 'groups': len(unique)}
    if len(unique) < 2:
        return {**result, 'lo': None, 'hi': None, 'valid_replicates': 0}
    sums = np.bincount(inverse, weights=values)
    counts = np.bincount(inverse)
    rng = np.random.default_rng(seed)
    draws = rng.integers(len(unique), size=(repetitions, len(unique)))
    estimates = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    lo, hi = np.quantile(estimates, [.025, .975])
    return {**result, 'lo': float(lo), 'hi': float(hi),
            'valid_replicates': len(estimates)}


def distance_bin(value):
    if not np.isfinite(value):
        return 'missing'
    if value < 0:
        raise ValueError('Negative genetic distance')
    if value <= 1e-4:
        return '<=0.0001'
    if value <= 1e-3:
        return '(0.0001,0.001]'
    if value <= 1e-2:
        return '(0.001,0.01]'
    return '>0.01'
