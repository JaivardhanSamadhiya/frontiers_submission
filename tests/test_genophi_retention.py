"""Test the exact upstream grid with mocked fits, not a reimplemented selector."""
import ast
import itertools
import logging
import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.benchmarks.genophi_retention import build_retaining_grid_search


def upstream_source():
    text = (ROOT/'data/external_validation/genophi_comparison/official_source/genophi/feature_selection.py').read_text(encoding='utf-8')
    func = next(n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == 'grid_search')
    return ast.get_source_segment(text, func)


def exercise(tmp_path, scores, repaired):
    seen = []
    class Model:
        def __init__(self, i): self.i = i
        def predict_proba(self, X): return np.array([[.8,.2],[.2,.8]])
    def fit(*args, **kwargs):
        i = args[5]['candidate']; seen.append(i)
        return Model(i), .5, .5, scores[i], np.array([0,1])
    noop = lambda *a, **k: None
    scope = dict(os=os, itertools=itertools, pd=pd, logging=logging,
                 train_and_evaluate=fit, plot_confusion_matrix=noop,
                 plot_roc_curve=noop, plot_precision_recall_curve=noop)
    source = upstream_source()
    if repaired: function = build_retaining_grid_search(source, scope)
    else:
        exec(compile(source, '<original_upstream_grid>', 'exec'), scope)
        function = scope['grid_search']
    X = pd.DataFrame({'feature':[0,1]}); y = pd.Series([0,1])
    ids = pd.DataFrame({'strain':['a','b'], 'phage':['p','p']})
    result = function(X,y,X,y,ids,ids,{'candidate':list(range(len(scores)))},str(tmp_path))
    assert seen == list(range(len(scores)))
    return result


@pytest.mark.parametrize('scores', [[-.7,-.1,-.4], [-.2,-.1,-.1]])
def test_negative_candidates_retained_with_original_last_tie_rule(tmp_path, scores):
    old = exercise(tmp_path/'old', scores, False)
    new = exercise(tmp_path/'new', scores, True)
    assert old[0] is None
    expected = max(i for i,s in enumerate(scores) if s == max(scores))
    assert new[0].i == expected and new[2] == max(scores)


@pytest.mark.parametrize('scores', [[-.4,0,-.1], [.4,.2,.4], [-.8,.1,.7]])
def test_nonnegative_maxima_preserve_original_result(tmp_path, scores):
    old = exercise(tmp_path/'old', scores, False)
    new = exercise(tmp_path/'new', scores, True)
    assert old[0].i == new[0].i and old[1:] == new[1:]


def test_unexpected_source_is_rejected():
    with pytest.raises(ValueError):
        build_retaining_grid_search('def grid_search():\n    best_mcc = 1\n', {})
