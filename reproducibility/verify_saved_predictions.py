"""Recalculate saved metrics and BASEL model scores without raw genomes."""
from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
ROOT = Path(__file__).resolve().parents[1]

def metrics(y, p):
    return [float(roc_auc_score(y,p)),float(average_precision_score(y,p))]

def main():
    result = {}
    folder = ROOT/'data/results_v2/basel_ecor_validation'
    scored = pd.read_csv(folder/'evaluated_predictions.csv')
    frozen = pd.read_csv(folder/'frozen_predictions.csv')
    matrix = np.load(folder/'frozen_input_matrix.npz',allow_pickle=False)['X']
    assert matrix.shape == (1240,444)
    model = joblib.load(ROOT/'data/results_v2/external_validation_models/sequence_only.joblib')
    recomputed = model['model'].predict_proba(model['scaler'].transform(matrix))[:,1]
    np.testing.assert_allclose(recomputed,frozen.score,atol=1e-7,rtol=0)
    joined = frozen.merge(scored,on=['host','phage'],validate='one_to_one',suffixes=('_frozen',''))
    np.testing.assert_allclose(joined.score_frozen,joined.score,atol=0,rtol=0)
    report = json.loads((folder/'evaluation.json').read_text())
    actual = metrics(joined.label,joined.score)
    np.testing.assert_allclose(actual,[report['metrics']['roc_auc'],report['metrics']['pr_auc']],atol=1e-13,rtol=0)
    result['basel_ecor'] = {'pairs':len(joined),'auroc':actual[0],'average_precision':actual[1],
                            'saved_model_scores_reproduced':True}
    staph = pd.read_csv(ROOT/'data/results_v2/external_staph_predictions.csv')
    saved = json.loads((ROOT/'data/results_v2/external_staph_validation.json').read_text())
    for name in ['full','sequence_only','vhip_four_feature']:
        actual = metrics(staph.label,staph['prob_'+name])
        expected = saved['models'][name]['pooled']
        np.testing.assert_allclose(actual,[expected['roc_auc'],expected['pr_auc']],atol=1e-13,rtol=0)
        result['staph_'+name] = {'auroc':actual[0],'average_precision':actual[1]}
    canonical = pd.read_csv(ROOT/'data/results_v2/canonical_sensitivity/predictions.csv')
    summaries = pd.read_csv(ROOT/'data/results_v2/canonical_sensitivity/summary.csv')
    for row in summaries.itertuples():
        part=canonical[(canonical.regime==row.regime)&(canonical.seed==row.seed)]
        assert len(part)==398 and part.row_index.is_unique
        np.testing.assert_allclose(metrics(part.label,part.score),[row.auroc,row.average_precision],atol=1e-13,rtol=0)
    result['canonical_conditions_recomputed'] = len(summaries)
    print(json.dumps(result,indent=2))
    print('PASS: saved-score/model checks only; raw-genome reconstruction not performed')

if __name__=='__main__':main()
