"""One-shot pipeline tail: wait for its fitting process, then independently audit.

Not a recurring monitor; does not push, submit, or modify manuscript claims.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json
import subprocess
import sys
import psutil

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT/'data/results_v2/genophi_comparison/official_modeling_full_budget_v1'


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--wait-for-pid', type=int)
    mode.add_argument('--run', action='store_true', help='Run fitting from verified checkpoints, then audit')
    args = parser.parse_args()
    try:
        if args.run:
            for process in psutil.process_iter(['pid', 'name', 'cmdline']):
                if 'python' in (process.info['name'] or '').lower() and any(
                        '47_official_genophi_outer_comparison.py' in c for c in (process.info['cmdline'] or [])):
                    raise RuntimeError('A fitting process is already active; do not start a duplicate')
            (RUN/'completion_status.json').write_text(json.dumps({'state': 'fitting_then_audit_pending',
                'updated_utc': datetime.now(timezone.utc).isoformat()}, indent=2))
            subprocess.run([sys.executable, str(ROOT/'experiments/47_official_genophi_outer_comparison.py')], cwd=ROOT, check=True)
        else:
            try:
                process = psutil.Process(args.wait_for_pid)
                if not any('47_official_genophi_outer_comparison.py' in c for c in process.cmdline()):
                    raise ValueError('Not the expected fitting process')
                process.wait()
            except psutil.NoSuchProcess:
                pass
        state = json.loads((RUN/'status.json').read_text())
        if state['state'] != 'requested_folds_predicted_pending_independent_audit' or state['folds'] != list(range(5)):
            raise RuntimeError('Fitting did not complete all five folds; audit not started')
        with (RUN/'independent_audit.log').open('w') as log:
            subprocess.run([sys.executable, str(ROOT/'experiments/49_audit_official_genophi_outer.py')], cwd=ROOT, stdout=log, stderr=log, check=True)
        audit = json.loads((RUN/'audit.json').read_text())
        if not audit['passed'] or audit['selected_models_reloaded'] != 250: raise ValueError('Audit not complete')
        state = {'state': 'complete_and_audited', 'report': str(ROOT/'GENOPHI_OFFICIAL_FULL_BUDGET_RESULTS.md')}
    except Exception as e:
        state = {'state': 'incomplete_or_failed', 'error': repr(e)}
    state['updated_utc'] = datetime.now(timezone.utc).isoformat()
    (RUN/'completion_status.json').write_text(json.dumps(state, indent=2))
    print(json.dumps(state), flush=True)
    if state['state'] != 'complete_and_audited': raise SystemExit(1)


if __name__ == '__main__': main()
