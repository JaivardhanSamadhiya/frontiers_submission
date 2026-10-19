"""Disclosed GenoPHI repair: retain a candidate even if every MCC is negative.

The pinned upstream files remain unchanged. Only the initial best_mcc constant
in the classification grid-search function is replaced in memory. See the
dated amendment; this is an adapted comparator, not unmodified GenoPHI.
"""
import ast
import hashlib
import inspect
import json
from pathlib import Path


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_retaining_grid_search(source, namespace):
    tree = ast.parse(source)
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    if len(functions) != 1 or functions[0].name != 'grid_search':
        raise ValueError('Expected exactly the upstream grid_search function')
    initializers = [n for n in functions[0].body if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == 'best_mcc' for t in n.targets)]
    if len(initializers) != 1 or not isinstance(initializers[0].value, ast.Constant) or initializers[0].value.value != 0:
        raise ValueError('Unexpected upstream best_mcc initialization')
    initializers[0].value = ast.parse("-float('inf')", mode='eval').body
    ast.fix_missing_locations(tree)
    local = {}
    exec(compile(tree, '<disclosed_genophi_negative_mcc_retention>', 'exec'), namespace, local)
    return local['grid_search']


def install_retention_repair():
    import genophi.feature_selection as selection
    import genophi.select_feature_modeling as modeling
    original = selection.grid_search
    repaired = build_retaining_grid_search(inspect.getsource(original), original.__globals__)
    selection.grid_search = repaired
    modeling.grid_search = repaired


def validate_repair_manifest(root, run):
    record = json.loads((run/'negative_mcc_repair_manifest.json').read_text())
    for name, expected in record['current_files_sha256'].items():
        if sha(root/name) != expected: raise ValueError(f'Repair source changed: {name}')
    for name, expected in record['preserved_files_sha256'].items():
        if sha(root/name) != expected: raise ValueError(f'Preserved evidence changed: {name}')
    earlier = json.loads((run/'execution_amendment.json').read_text())
    manifest = json.loads((run/'run_manifest.json').read_text())
    driver_key = str(Path('experiments')/'47_official_genophi_outer_comparison.py')
    if manifest['sha256'][driver_key] != earlier['old_driver_sha256']:
        raise ValueError('Original driver identity changed')
    if sha(run/'driver_before_prediction_api_fix.py') != earlier['old_driver_sha256']:
        raise ValueError('Original driver snapshot changed')
    if sha(run/'before_negative_mcc_repair/47_official_genophi_outer_comparison.py') != earlier['new_driver_sha256']:
        raise ValueError('Repair chain broken')
    return record


def verify_model_selection(directory):
    """Verify retained model matches the final maximizer, including MCC ties."""
    import joblib
    import numpy as np
    import pandas as pd
    from sklearn.metrics import matthews_corrcoef
    directory = Path(directory)
    table = pd.read_csv(directory/'model_performance.csv')
    if len(table) != 8 or not np.isfinite(table.mcc).all():
        raise ValueError('Incomplete/nonfinite parameter search')
    best = table[table.mcc == table.mcc.max()].iloc[-1]
    model = joblib.load(directory/'best_model.pkl')
    params = model.get_params()
    for key in ['iterations', 'learning_rate', 'depth', 'loss_function', 'thread_count']:
        if params[key] != best[key]: raise ValueError(f'Model does not match grid choice: {key}')
    predictions = pd.read_csv(directory/'best_model_predictions.csv')
    mcc = matthews_corrcoef(predictions.interaction, predictions.Prediction)
    if not np.isclose(mcc, best.mcc, rtol=0, atol=1e-12):
        raise ValueError('Selected predictions do not match selected MCC')
    return {'grid_candidates': 8, 'best_development_mcc': float(best.mcc),
            'all_candidates_negative': bool((table.mcc < 0).all()),
            'selection': 'Last maximum MCC in official grid order; no sign-based rejection',
            'model_sha256': sha(directory/'best_model.pkl'),
            'performance_sha256': sha(directory/'model_performance.csv'),
            'predictions_sha256': sha(directory/'best_model_predictions.csv')}
