"""Batch-invariant ranking of a fixed bank using development reference scores.

The output is a ranking score, not a calibrated probability or biological
independence model. Fit reference distributions on development predictions only.
"""
from __future__ import annotations

import numpy as np


class DevelopmentRankNormalizer:
    def fit(self, reference):
        reference = np.asarray(reference, dtype=float)
        if reference.ndim != 2 or not reference.shape[0] or not reference.shape[1]:
            raise ValueError("Reference must be a non-empty hosts-by-phages matrix")
        if not np.isfinite(reference).all():
            raise ValueError("Reference scores must be finite")
        self.reference_ = np.sort(reference, axis=0)
        return self

    def transform(self, scores):
        if not hasattr(self, "reference_"):
            raise RuntimeError("Fit development references before scoring queries")
        scores = np.asarray(scores, dtype=float)
        if scores.ndim != 2 or scores.shape[1] != self.reference_.shape[1]:
            raise ValueError("Query columns must match the fitted phage order")
        if not np.isfinite(scores).all():
            raise ValueError("Query scores must be finite")
        result = np.empty_like(scores)
        n = len(self.reference_)
        for j in range(scores.shape[1]):
            reference = self.reference_[:, j]
            left = np.searchsorted(reference, scores[:, j], side="left")
            right = np.searchsorted(reference, scores[:, j], side="right")
            result[:, j] = (left + right + 1) / (2 * (n + 1))
        return result


def combined_ranking_score(constituent_scores):
    scores = np.asarray(constituent_scores, dtype=float)
    if scores.ndim != 2 or not scores.shape[1]:
        raise ValueError("Supply a hosts-by-constituents matrix")
    if not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError("Constituent ranking scores must lie in [0, 1]")
    return 1 - np.prod(1 - scores, axis=1)
