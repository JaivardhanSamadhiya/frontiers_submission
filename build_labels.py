"""
build_labels.py - Construct negative phage-host pairs and validate splits
==========================================================================
- DataLeakageError: raised when a pair appears in both train and test sets
- build_negatives(pos_df, cfg, seed): generates within-genus + cross-genus
  negatives and returns the full positive+negative dataset.
- validate_no_leakage(train_df, test_df, context): per-fold sanity check.

This module is local-only: no network calls, no external state.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class DataLeakageError(Exception):
    """Raised when train and test splits share a phage-host pair."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def build_negatives(pos_df: pd.DataFrame, cfg: dict, seed: int) -> pd.DataFrame:
    """Generate within-genus and cross-genus negatives from positives only.

    Returns the combined dataset (positives + negatives).  All negatives carry
    is_constructed_negative = True; all positives carry False.

    Numeric columns that exist in pos_df but not on a generated negative are
    imputed using the 75th percentile of the positive column.  Missing columns
    are imputed with 0.0 (with a warning).
    """
    rng = np.random.default_rng(seed)

    pos = pos_df.copy()
    if "label" not in pos.columns:
        pos["label"] = 1
    if "is_constructed_negative" not in pos.columns:
        pos["is_constructed_negative"] = False
    if "host_genus" not in pos.columns:
        pos["host_genus"] = pos["host"].astype(str).str.split(" ").str[0]
    if "source" not in pos.columns:
        pos["source"] = "unknown"

    # Index positive pairs for quick lookup
    pos_pairs = set(zip(pos["phage"], pos["host"]))

    all_hosts = sorted(pos["host"].unique().tolist())
    host_to_genus = dict(zip(pos["host"], pos["host_genus"]))
    genera = sorted(set(host_to_genus.values()))

    # for each genus, list of hosts in it
    genus_to_hosts: dict[str, list[str]] = {}
    for h in all_hosts:
        g = host_to_genus.get(h, "")
        genus_to_hosts.setdefault(g, []).append(h)

    # for each phage, hosts and genera it positively infects
    phage_pos_hosts: dict[str, set[str]] = {}
    phage_pos_genera: dict[str, set[str]] = {}
    for _, r in pos.iterrows():
        phage_pos_hosts.setdefault(r["phage"], set()).add(r["host"])
        phage_pos_genera.setdefault(r["phage"], set()).add(host_to_genus.get(r["host"], ""))

    n_within = int(cfg["data"]["neg_ratio_within_genus"])
    n_cross = int(cfg["data"]["neg_ratio_cross_genus"])

    new_rows: list[dict] = []

    for phage, pos_hosts in phage_pos_hosts.items():
        pos_gen = phage_pos_genera[phage]

        # Within-genus negatives
        candidates = []
        for g in pos_gen:
            candidates.extend([h for h in genus_to_hosts.get(g, []) if h not in pos_hosts])
        candidates = list(set(candidates))
        if candidates:
            k = min(n_within * len(pos_hosts), len(candidates))
            if k > 0:
                chosen = rng.choice(candidates, size=k, replace=False)
                for h in chosen:
                    pair = (phage, h)
                    if pair in pos_pairs:
                        continue
                    new_rows.append({
                        "phage": phage,
                        "host": h,
                        "host_genus": host_to_genus.get(h, ""),
                        "label": 0,
                        "is_constructed_negative": True,
                        "source": "constructed_within_genus",
                    })

        # Cross-genus negatives
        non_pos_genera = [g for g in genera if g not in pos_gen and g]
        if non_pos_genera:
            cand_hosts = []
            for g in non_pos_genera:
                cand_hosts.extend(genus_to_hosts.get(g, []))
            cand_hosts = list(set(cand_hosts) - pos_hosts)
            k = min(n_cross * len(pos_hosts), len(cand_hosts))
            if k > 0:
                chosen = rng.choice(cand_hosts, size=k, replace=False)
                for h in chosen:
                    pair = (phage, h)
                    if pair in pos_pairs:
                        continue
                    new_rows.append({
                        "phage": phage,
                        "host": h,
                        "host_genus": host_to_genus.get(h, ""),
                        "label": 0,
                        "is_constructed_negative": True,
                        "source": "constructed_cross_genus",
                    })

    neg_df = pd.DataFrame(new_rows)
    log.info(f"[negatives] generated {len(neg_df)} candidate negatives "
             f"(within={int((neg_df.get('source') == 'constructed_within_genus').sum()) if len(neg_df) else 0}, "
             f"cross={int((neg_df.get('source') == 'constructed_cross_genus').sum()) if len(neg_df) else 0})")

    # Numeric imputation for negatives - use 75th percentile of positives
    numeric_cols = pos.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [c for c in numeric_cols if c not in
                    ("label",)]
    if len(neg_df) > 0:
        for col in numeric_cols:
            if col in pos.columns:
                try:
                    q75 = float(np.nanpercentile(pos[col].astype(float).values, 75))
                except Exception:
                    q75 = 0.0
                    log.warning(f"[negatives] could not compute q75 for {col}; imputing 0.0")
            else:
                q75 = 0.0
                log.warning(f"[negatives] column '{col}' missing in positives; imputing 0.0")
            neg_df[col] = q75

    # Combine
    out_cols = list(set(list(pos.columns) + list(neg_df.columns)))
    pos_out = pos.reindex(columns=out_cols)
    neg_out = neg_df.reindex(columns=out_cols)
    combined = pd.concat([pos_out, neg_out], ignore_index=True)

    # Drop duplicates (positives win over negatives)
    combined = combined.sort_values("label", ascending=False).drop_duplicates(
        ["phage", "host"], keep="first").reset_index(drop=True)

    log.info(f"[negatives] final dataset: {len(combined)} rows "
             f"(pos={int((combined['label'] == 1).sum())}, "
             f"neg={int((combined['label'] == 0).sum())})")

    # Save a methodology report
    try:
        results_dir = Path(cfg["paths"]["results_dir"])
    except Exception:
        results_dir = Path("data/results")
    results_dir.mkdir(parents=True, exist_ok=True)
    report = results_dir / "negative_sampling_report.txt"
    report.write_text(
        "Negative pair construction\n"
        "==========================\n"
        f"Seed: {seed}\n"
        f"Within-genus ratio: {n_within} (per positive host)\n"
        f"Cross-genus ratio:  {n_cross} (per positive host)\n"
        f"Total positives: {int((combined['label'] == 1).sum())}\n"
        f"Total negatives: {int((combined['label'] == 0).sum())}\n"
        "\nMethodology\n"
        "-----------\n"
        "For each phage p with observed positive hosts P_p:\n"
        "  1. Within-genus negatives: sample hosts from the same bacterial\n"
        "     genus as members of P_p that are NOT in P_p.\n"
        "  2. Cross-genus negatives: sample hosts from genera not associated\n"
        "     with any positive host for p.\n"
        "Constructed negatives are flagged is_constructed_negative=True.\n"
        "Numeric feature values on negatives are imputed using the 75th\n"
        "percentile of the corresponding positive distribution.\n"
        "\nCaveat: absence of recorded interaction does not guarantee\n"
        "biological non-infectivity (open-world assumption).\n",
        encoding="utf-8",
    )
    log.info(f"[negatives] methodology report -> {report}")
    return combined


def validate_no_leakage(train_df: pd.DataFrame, test_df: pd.DataFrame,
                        context: str = "") -> None:
    """Raise DataLeakageError if any (phage, host) pair appears in both."""
    train_pairs = set(zip(train_df["phage"].astype(str).str.lower().str.strip(),
                          train_df["host"].astype(str).str.lower().str.strip()))
    test_pairs = set(zip(test_df["phage"].astype(str).str.lower().str.strip(),
                         test_df["host"].astype(str).str.lower().str.strip()))
    overlap = train_pairs & test_pairs
    if overlap:
        raise DataLeakageError(
            f"[{context}] {len(overlap)} pairs appear in both train and test. "
            f"First 5: {list(overlap)[:5]}")
    log.debug(f"[{context}] Leakage check passed.")
