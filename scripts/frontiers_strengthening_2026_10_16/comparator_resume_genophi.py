"""Resume the historical, disclosed GenoPHI adapter in an isolated run copy.

No modeling code, budgets, inputs, or historical checkpoint is changed. Existing
checkpoint hashes continue to verify their preserved original files. Newly
completed stages point to this run directory. Subprocess calls are redirected to
this launcher so workers also use the isolated run directory.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / 'data/results_v2/genophi_comparison/official_modeling_full_budget_v1'
OUT = ROOT / 'data/results_v2/frontiers_strengthening_2026_10_16/comparator'
RUN = OUT / 'official_full_budget_resumed'
DRIVER = ROOT / 'experiments/47_official_genophi_outer_comparison.py'
AUDITOR = ROOT / 'experiments/49_audit_official_genophi_outer.py'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.RUN = RUN
    return module


def prepare_copy():
    if (OUT / 'isolated_resume_provenance.json').exists():
        record = json.loads((OUT / 'isolated_resume_provenance.json').read_text())
        assert digest(DRIVER) == record['driver_sha256']
        return
    OUT.mkdir(parents=True, exist_ok=True)
    original_hashes = {}
    for count, path in enumerate(p for p in ORIGINAL.rglob('*') if p.is_file()):
        name = str(path.relative_to(ORIGINAL))
        value = digest(path)
        original_hashes[name] = value
        destination = RUN / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Resume a partial initial snapshot without recopying matching files.
        if not destination.exists() or digest(destination) != value:
            shutil.copy2(path, destination)
        assert digest(destination) == value, name
        if count % 1000 == 0:
            print(json.dumps({'stage': 'verify_isolated_copy', 'files': count + 1}), flush=True)
    record = {
        'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'source': str(ORIGINAL), 'isolated_destination': str(RUN),
        'source_files_sha256': original_hashes,
        'driver_sha256': digest(DRIVER), 'auditor_sha256': digest(AUDITOR),
        'launcher_sha256': digest(Path(__file__)),
        'scope': 'Full-budget adapted GenoPHI on previously used Pseudomonas outer folds; not a new untouched cohort or exact unmodified reproduction.',
        'unchanged': ['25 selection runs', '50 modeling runs per eligible cutoff',
                      'official 8-setting grid', '4 threads', 'negative-MCC retention repair',
                      'outer folds', 'source-only protein reference construction'],
    }
    (OUT / 'isolated_resume_provenance.json').write_text(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', choices=['prepare', 'selection', 'tables', 'model', 'predict'])
    parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--repeat', type=int)
    parser.add_argument('--cutoff', type=int)
    parser.add_argument('--audit', action='store_true')
    args = parser.parse_args()
    if args.worker:
        load(DRIVER, 'isolated_genophi_driver').worker(args)
        return
    prepare_copy()
    if args.audit:
        module = load(AUDITOR, 'isolated_genophi_auditor')
        # Only report destination changes; numerical auditor source is retained.
        source = AUDITOR.read_text(encoding='utf-8')
        target = "(ROOT/'GENOPHI_OFFICIAL_FULL_BUDGET_RESULTS.md')"
        assert source.count(target) == 1
        source = source.replace(target, "(RUN/'GENOPHI_OFFICIAL_FULL_BUDGET_RESULTS.md')")
        scope = {'__name__': 'isolated_genophi_auditor', '__file__': str(AUDITOR)}
        exec(compile(source, str(AUDITOR), 'exec'), scope)
        scope['RUN'] = RUN
        scope['main']()
        return
    driver = load(DRIVER, 'isolated_genophi_driver')
    original_run = subprocess.run

    def redirected_run(command, *positional, **keyword):
        command = [str(Path(__file__)) if str(item) == str(DRIVER) else item for item in command]
        return original_run(command, *positional, **keyword)

    subprocess.run = redirected_run
    driver.main(args)
    subprocess.run = original_run
    original_run([sys.executable, str(Path(__file__)), '--audit'], cwd=ROOT, check=True)
    (OUT / 'completion.json').write_text(json.dumps({
        'state': 'all_five_folds_complete_and_audited',
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'audit': str(RUN / 'audit.json'),
    }, indent=2))


if __name__ == '__main__':
    main()
