"""Read-only feasibility audit for the *official* Gaborieau comparator.

Does not read external cocktail outcomes, fit models, fill missing scores,
unpickle downloaded objects, or modify upstream source files. Outputs are audit
records only. The official training script is dry-run under intercepted mkdir
and read_csv so that the first missing hard-coded input is documented safely.
"""
from __future__ import annotations

import ast
import concurrent.futures
import contextlib
import hashlib
import io
import json
import platform
import re
import subprocess
import traceback
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from unittest.mock import patch
from urllib.request import Request, urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "external/gaborieau_coli_phage_interactions_2023"
OUT = ROOT / "data/results_v2/frontiers_revision_2026_10_15/official_gaborieau"
PIN = "99656d8c020fcd1dbbc30a02ee10274a07655562"
API = "https://api.github.com/repos/mdmparis/coli_phage_interactions_2023"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(
        ["git", "-c", f"safe.directory={SOURCE.as_posix()}", "-C", str(SOURCE), *args],
        text=True,
    ).strip()


def get_json(url):
    request = Request(url, headers={"User-Agent": "PrecisionPhage-Frontiers-revision-audit/1.0"})
    with urlopen(request, timeout=60) as result:
        return json.load(result)


def archive_member_hosts(member):
    # RAR is never extracted to filesystem. tar streams one named, validated
    # ordinary CSV member; only identifiers/shape are retained in the report.
    archive = SOURCE / "dev/predictions/results/predictions.rar"
    result = subprocess.run(["tar", "-xOf", str(archive), member],
                            capture_output=True, check=True)
    frame = pd.read_csv(io.BytesIO(result.stdout), sep=";", usecols=["bacteria", "dataset"])
    return {"member": member, "rows": len(frame),
            "hosts": sorted(frame.bacteria.unique().tolist()),
            "datasets": sorted(frame.dataset.unique().tolist()),
            "sha256": hashlib.sha256(result.stdout).hexdigest()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert (OUT / "PROTOCOL.md").exists(), "Write the protocol before running this audit"
    files = [p for p in SOURCE.rglob("*") if p.is_file() and ".git" not in p.parts]
    original = {p.relative_to(SOURCE).as_posix(): sha(p) for p in files}
    assert git("rev-parse", "HEAD") == PIN
    assert not git("status", "--short"), "Upstream checkout is not pristine"
    current = get_json(API + "/commits/main")
    tree = get_json(API + "/git/trees/" + current["sha"] + "?recursive=1")
    (OUT / "upstream_tree.json").write_text(json.dumps(tree, indent=2), encoding="utf-8")
    (OUT / "upstream_commit.json").write_text(json.dumps(current, indent=2), encoding="utf-8")
    remote_files = [x["path"] for x in tree["tree"] if x["type"] == "blob"]
    executable = [x for x in remote_files if x.endswith((".py", ".pkl", ".pickle", ".joblib", ".onnx"))]

    step = pd.read_csv(SOURCE / "dev/cocktails/data/cocktails_step.csv", sep=";")
    target = pd.read_csv(SOURCE / "dev/cocktails/data/test_collection.csv", sep=";")
    train = pd.read_csv(SOURCE / "data/genomics/bacteria/picard_collection.csv", sep=";")
    umap = pd.read_csv(SOURCE / "data/genomics/bacteria/umap_phylogeny/coli_umap_8_dims.tsv", sep="\t")
    selected_hosts = set(target.bacteria)
    missing = step.loc[step.y_pred_proba.isna(), ["bacteria", "phage", "chosen_at_step", "rank"]].copy()
    missing.to_csv(OUT / "missing_official_constituent_score_request.csv", index=False)
    stage = step.groupby("chosen_at_step").agg(rows=("phage", "size"), numeric_scores=("y_pred_proba", "count"))
    stage.to_csv(OUT / "score_coverage_by_stage.csv")
    assert len(step) == 300 and len(missing) == 67
    assert set(missing.chosen_at_step) == {"4_rescue"}
    assert not step[["bacteria", "phage"]].duplicated().any()
    complete = step.groupby("bacteria").y_pred_proba.count()
    assert (complete == 3).sum() == 64

    # Inspect the official script's source and safely observe its first failed
    # load without creating the authors' hard-coded directories or training.
    official = SOURCE / "dev/predictions/predict_all_phages.py"
    code = official.read_text(encoding="utf-8")
    hardcoded_inputs = []
    for node in ast.walk(ast.parse(code)):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "pd"
                and node.func.attr == "read_csv" and node.args
                and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str)):
            literal_path = node.args[0].value
            hardcoded_inputs.append({"line": node.lineno, "path": literal_path,
                                     "exists_in_current_environment": Path(literal_path).is_file()})
    namespace = {"__name__": "__official_preflight__", "__file__": str(official)}
    requested_reads, requested_mkdirs = [], []

    def no_mkdir(path, *args, **kwargs):
        requested_mkdirs.append(str(path))

    def stop_read(path, *args, **kwargs):
        requested_reads.append(str(path))
        assert not Path(path).is_file(), "Unexpected private source file exists; stop and inspect"
        raise FileNotFoundError(f"Official source input does not exist: {path}")

    stream = io.StringIO()
    with contextlib.redirect_stdout(stream), patch("os.mkdir", no_mkdir), patch("pandas.read_csv", stop_read):
        try:
            exec(compile(code, str(official), "exec"), namespace)
        except FileNotFoundError as exc:
            dry_run = {"exception": type(exc).__name__, "message": str(exc),
                       "requested_reads": requested_reads, "intercepted_output_mkdir": requested_mkdirs,
                       "stdout": stream.getvalue(), "fitting_started": False,
                       "source_executed_unmodified": True,
                       "side_effect_interception": "os.mkdir suppressed; pandas.read_csv stopped before read"}
        else:
            raise AssertionError("Expected explicit missing-input stop did not occur")

    # Independently scan every released prediction CSV, not merely one example.
    archive = SOURCE / "dev/predictions/results/predictions.rar"
    listed = subprocess.check_output(["tar", "-tf", str(archive)], text=True).splitlines()
    assert listed and all(not PurePosixPath(x).is_absolute() and ".." not in PurePosixPath(x).parts
                          for x in listed)
    assert all(x == "predictions" or x.endswith(".csv") for x in listed)
    names = [x for x in listed if x.endswith(".csv")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        members = list(pool.map(archive_member_hosts, names))
    archive_hosts = sorted(set().union(*(set(x["hosts"]) for x in members)))
    target_matches = sorted(selected_hosts & set(archive_hosts))
    for member in members:
        member["target_name_overlap"] = sorted(selected_hosts & set(member.pop("hosts")))
    (OUT / "prediction_archive_inventory.json").write_text(
        json.dumps({"members": members, "unique_hosts": archive_hosts,
                    "target_name_matches": target_matches,
                    "warning": "AN03 is a known development/query identifier collision, not an external prediction"}, indent=2),
        encoding="utf-8")

    # No remote issue or email is sent; existing public issues can be checked.
    issues = get_json(API + "/issues?state=all&per_page=100")
    issue_summary = [{"number": x["number"], "title": x["title"], "url": x["html_url"],
                      "state": x["state"], "body": x.get("body"), "comments": x["comments"],
                      "comment_details": get_json(x["comments_url"]) if x["comments"] else []}
                     for x in issues]
    (OUT / "upstream_existing_issues.json").write_text(json.dumps(issue_summary, indent=2), encoding="utf-8")
    unchanged = all(sha(SOURCE / p) == expected for p, expected in original.items())
    assert unchanged and not git("status", "--short")

    report = {
        "status": "faithful all-100 official predictor rerun not executable from released resources",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": sha(OUT / "PROTOCOL.md"),
        "source_commit_pinned": PIN, "source_current_main": current["sha"],
        "source_current_equals_pinned": current["sha"] == PIN,
        "remote_tree_entries": len(tree["tree"]), "upstream_file_count": len(files),
        "remote_executable_or_model_files": executable,
        "fitted_classifier_files_found": [],
        "only_pickle_is_unrelated_defense_data_not_executed": True,
        "upstream_source_sha256": original,
        "score_coverage": {"selected_constituents": len(step), "numeric": int(step.y_pred_proba.count()),
                           "missing": len(missing), "complete_hosts": int((complete == 3).sum()),
                           "hosts_with_missing": int((complete < 3).sum()),
                           "missing_stage": "4_rescue", "stages": stage.reset_index().to_dict("records")},
        "external_input_coverage": {"hosts": len(target), "target_feature_columns": list(target.columns),
                                     "missing_required_capsule_columns": ["Klebs_capsule_type", "ABC_serotype"],
                                     "published_training_umap_hosts": len(umap),
                                     "target_names_overlapping_published_umap": sorted(selected_hosts & set(umap.bacteria)),
                                     "target_names_overlapping_picard": sorted(selected_hosts & set(train.bacteria)),
                                     "AN03_is_identifier_collision_not_valid_external_feature_match": True},
        "official_script_safe_preflight": dry_run,
        "official_script_hardcoded_input_paths": hardcoded_inputs,
        "prediction_archive": {"members_checked": len(members), "total_rows": sum(x["rows"] for x in members),
                               "unique_development_host_names": len(archive_hosts),
                               "target_name_overlap": target_matches,
                               "authentic_external_host_predictions_recovered": 0},
        "blocking_resources": [
            "No released fitted pointwise classifier, preprocessing object, or executable four-stage cocktail recommender.",
            "The official pointwise CV script references missing private Windows inputs and the README explicitly says the original CV cluster file is unavailable.",
            "The complete external capsule feature table and external eight-dimensional UMAP transform/coordinates are not released.",
            "All 67 unavailable values are rescue selections. The source README defines rescue by generalist coverage, not a numeric host-specific predicted probability.",
            "The released prediction RAR contains development cross-validation outputs, not the 100 external-host outputs.",
        ],
        "interpretation": "Do not fill missing rescue values or call a newly reconstructed classifier the official frozen comparator. Keep the exact 64-complete-host common comparison and ask authors whether an all-100 continuous score exists.",
        "new_outcome_evaluation_performed": False, "external_cocktail_outcome_files_read": False,
        "models_fitted": 0, "gpu_used": False, "emails_or_issues_sent": False,
        "source_files_unchanged": unchanged,
        "software": {"python": platform.python_version(), "pandas": pd.__version__},
        "public_sources": {"repository": "https://github.com/mdmparis/coli_phage_interactions_2023/tree/" + PIN,
                           "article": "https://doi.org/10.1038/s41564-024-01832-5",
                           "contact_verification": "https://pubmed.ncbi.nlm.nih.gov/39482383/"},
    }
    (OUT / "official_comparator_feasibility.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in ["status", "score_coverage", "external_input_coverage",
                                            "prediction_archive", "source_files_unchanged"]}, indent=2))


if __name__ == "__main__":
    main()
