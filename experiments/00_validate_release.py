#!/usr/bin/env python3
"""Fast, read-only validation of frozen inputs and submission artifacts."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from precisionphage.data.genomes import GenomeIndex  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix.lower() in {".csv", ".json"}:
        # Git may check text artifacts out with CRLF or LF. Hash a canonical LF
        # representation so the release manifest is platform-independent.
        text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
        digest.update(text.encode("utf-8"))
        return digest.hexdigest()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-missing-genomes", action="store_true",
        help="permit a public code-only clone without the large local FASTA bundle")
    args = parser.parse_args()
    raw = pd.read_csv(ROOT / "data/raw/VirusHostInter.csv")
    modelable = pd.read_csv(ROOT / "data/interim_v2/interactions_modelable.csv")
    require(len(raw) == 8849, f"raw row count changed: {len(raw)}")
    require(len(modelable) == 1947,
            f"modelable row count changed: {len(modelable)}")
    require(int(modelable["label"].sum()) == 1488,
            "modelable positive count changed")
    require(set(modelable["study"]) == {"NCBI_HR"},
            "sequence-covered data are no longer the documented NCBI_HR subset")

    phage_dir = ROOT / "external/phist_run/phages"
    host_dir = ROOT / "external/phist_run/hosts"
    genomes_present = (phage_dir.is_dir() and host_dir.is_dir()
                       and any(phage_dir.glob("*.fasta"))
                       and any(host_dir.glob("*.fasta")))
    if genomes_present:
        phage_index = GenomeIndex([phage_dir])
        host_index = GenomeIndex([host_dir])
        phage_cov = phage_index.coverage(sorted(modelable["phage"].unique()))
        host_cov = host_index.coverage(sorted(modelable["host"].unique()))
        require(phage_cov["resolved"] == phage_cov["n"],
                f"unresolved frozen phages: {phage_cov['missing'][:5]}")
        require(host_cov["resolved"] == host_cov["n"],
                f"unresolved frozen hosts: {host_cov['missing'][:5]}")
    else:
        require(args.allow_missing_genomes,
                "frozen FASTAs are absent; stage them or pass --allow-missing-genomes")
        print("NOTE: large frozen FASTAs absent; exact genome coverage check skipped")

    manifest_path = ROOT / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for rel, expected in manifest["sha256"].items():
        path = ROOT / rel
        require(path.is_file(), f"missing frozen artifact: {rel}")
        observed = sha256(path)
        require(observed == expected,
                f"checksum mismatch for {rel}: {observed} != {expected}")

    required = [
        "data/results_v2/leakage_splits_results.json",
        "data/results_v2/gnn_ablation.json",
        "data/results_v2/significance_results.json",
        "data/results_v2/cocktail_summary.json",
        "data/results_v2/temporal_summary.json",
        "data/results_v2/temporal_trajectory.npz",
    ]
    for rel in required:
        require((ROOT / rel).is_file(), f"missing required artifact: {rel}")
    coverage = "exact genome coverage" if genomes_present else "genome check skipped"
    print(f"PASS: frozen inputs, {coverage}, checksums, and required artifacts")


if __name__ == "__main__":
    main()
