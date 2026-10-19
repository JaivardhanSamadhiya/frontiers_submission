"""One-time verification/registration before resuming the adapted comparator."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.benchmarks.genophi_retention import sha, verify_model_selection
RUN = ROOT/'data/results_v2/genophi_comparison/official_modeling_full_budget_v1'


def main():
    output = RUN/'negative_mcc_repair_manifest.json'
    if output.exists(): raise ValueError('Repair registration already exists; do not replace provenance')
    legacy, preserved = [], {}
    for marker in sorted(RUN.glob('0.001_fold*/stages/*.json')):
        status = json.loads(marker.read_text())
        if status['exit_code'] != 0: raise ValueError('Invalid existing checkpoint')
        preserved[str(marker.relative_to(ROOT))] = sha(marker)
        for name, expected in status['sha256'].items():
            if sha(ROOT/name) != expected: raise ValueError(f'Checkpoint mismatch: {name}')
            preserved[name] = expected
        if marker.stem.startswith('model_'):
            _, cutoff, repeat = marker.stem.split('_')
            directory = marker.parent.parent/f'modeling/cutoff_{cutoff}/run_{repeat}'
            result = verify_model_selection(directory)
            if result['best_development_mcc'] < 0: raise ValueError('Cannot reuse unequal selection semantics')
            for name in ['model_performance.csv','best_model_predictions.csv','feature_importances.csv']:
                path = directory/name; preserved[str(path.relative_to(ROOT))] = sha(path)
            legacy.append({'checkpoint': str(marker.relative_to(ROOT)), **result})
        elif marker.stem.startswith('selection_'):
            repeat = marker.stem.split('_')[-1]
            directory = marker.parent.parent/f'selection/run_{repeat}'
            table = pd.read_csv(directory/'model_performance.csv')
            preds = pd.read_csv(directory/'best_model_predictions.csv')
            if len(table) != 8 or not np.isfinite(table.mcc).all() or table.mcc.max() < 0:
                raise ValueError('Legacy feature-selection grid cannot be reused')
            if not np.isclose(matthews_corrcoef(preds.interaction,preds.Prediction),table.mcc.max(),atol=1e-12,rtol=0):
                raise ValueError('Legacy feature-selection choice mismatch')
            for name in ['model_performance.csv','best_model_predictions.csv']:
                path = directory/name; preserved[str(path.relative_to(ROOT))] = sha(path)
    if len(legacy) != 82: raise ValueError('Unexpected inherited model count')
    backups = RUN/'before_negative_mcc_repair'
    for path in backups.rglob('*'):
        if path.is_file(): preserved[str(path.relative_to(ROOT))] = sha(path)
    current = ['experiments/47_official_genophi_outer_comparison.py',
               'experiments/49_audit_official_genophi_outer.py',
               'scripts/finish_official_genophi_comparison.py',
               'src/precisionphage/benchmarks/genophi_retention.py',
               'tests/test_genophi_retention.py',
               'scripts/register_negative_mcc_repair.py',
               'GENOPHI_NEGATIVE_MCC_AMENDMENT_2026_10_06.md']
    record = {'registered_utc': datetime.now(timezone.utc).isoformat(),
              'retained_completed_models': len(legacy), 'retained_selection_repetitions': 25,
              'legacy_selection_equivalence': 'Each existing maximizing MCC is nonnegative; changing initialization leaves the last maximizer unchanged.',
              'scope': 'Adapted GenoPHI; negative-MCC retention repair, not unmodified official software.',
              'legacy_models': legacy,
              'current_files_sha256': {str(Path(p)):sha(ROOT/p) for p in current},
              'preserved_files_sha256': preserved}
    output.write_text(json.dumps(record,indent=2))
    print(json.dumps({'registered_utc':record['registered_utc'],'legacy_models_verified':len(legacy),
                      'preserved_checksums':len(preserved),'path':str(output)},indent=2))


if __name__ == '__main__': main()
