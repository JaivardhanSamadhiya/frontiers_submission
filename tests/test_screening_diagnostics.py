import itertools
import numpy as np
import pytest
from precisionphage.benchmarks.screening import (
    screening_curve, random_curve, cluster_mean_ci, distance_bin,
)


def test_stable_ties_and_all_negative():
    a = screening_curve([1, 0, 0], [.5, .5, .5], ['b', 'a', 'c'])
    assert list(a['order']) == ['a', 'b', 'c']
    assert a['first_hit'] == 2
    assert a['normalized_hit_area'] == 2/3
    empty = screening_curve([0, 0], [0., 1.], ['a', 'b'])
    assert np.isnan(empty['recall']).all()
    assert not empty['hit'].any()


def test_exact_random_matches_enumeration():
    for n, p in [(4, 2), (4, 0), (4, 4), (5, 1)]:
        curves = []
        for positives in itertools.combinations(range(n), p):
            y = np.zeros(n, dtype=int); y[list(positives)] = 1
            curves.append(screening_curve(y, -np.arange(n), list(map(str, range(n)))))
        random = random_curve(n, p)
        for metric in ['hit', 'precision', 'normalized_hit_area']:
            np.testing.assert_allclose(np.mean([c[metric] for c in curves], axis=0), random[metric])
        if p:
            assert np.mean([c['first_hit'] for c in curves]) == pytest.approx(random['first_hit'])


def test_bootstrap_missing_single_group_and_pairing():
    r = cluster_mean_ci([1, 1, 1, np.nan], ['a', 'b', 'b', 'c'])
    assert r['mean'] == r['lo'] == r['hi'] == 1
    assert r['n'] == 3 and r['groups'] == 2
    assert cluster_mean_ci([0, 1], ['a', 'a'])['lo'] is None
    assert cluster_mean_ci([np.nan], ['a'])['mean'] is None
    assert distance_bin(1e-4) == '<=0.0001'
    assert distance_bin(1e-3) == '(0.0001,0.001]'
    assert distance_bin(np.nan) == 'missing'


def test_reject_unknown_labels():
    with pytest.raises(ValueError):
        screening_curve([-1, 1], [.1, .8], ['a', 'b'])
