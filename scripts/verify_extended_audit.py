"""Check completed extended analyses against saved predictions and source labels."""
from pathlib import Path
import json
import zipfile
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'data/results_v2/canonical_sensitivity'
    pred=pd.read_csv(out/'predictions.csv')
    folds=pd.read_csv(out/'folds.csv')
    summary=pd.read_csv(out/'summary.csv')
    assert len(folds)==60 and (folds.status=='completed').all()
    assert len(summary)==12
    original=pd.read_csv(ROOT/'data/results_v2/post_rejection_audit/matched_cohort_predictions.csv')
    assert set(pred.row_index)==set(original.row_index)
    for (seed,fold),group in folds.groupby(['seed','fold']):
        assert group.train_n.nunique()==1 and group.train_positive.nunique()==1
    for r in summary.to_dict('records'):
        p=pred[(pred.seed==r['seed'])&(pred.regime==r['regime'])]
        assert len(p)==398 and not p.row_index.duplicated().any()
        assert np.isclose(roc_auc_score(p.label,p.score),r['auroc'])
        assert np.isclose(average_precision_score(p.label,p.score),r['average_precision'])
    out=ROOT/'data/results_v2/additional_banks'
    records=pd.read_csv(out/'host_selections.csv',dtype={'grouping':str})
    report=json.loads((out/'summary.json').read_text())
    with zipfile.ZipFile(ROOT/'data/external_validation/genophi_2026/processed_data.zip') as z:
        labels=pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'),dtype={'strain':str,'phage':str})
    for (bank,grouping),frame in records.groupby(['bank','grouping']):
        lookup=labels[labels.dataset==bank].set_index(['strain','phage']).interaction
        assert not frame.host.duplicated().any()
        assert frame.groupby('group').fold.nunique().max()==1
        for row in frame.to_dict('records'):
            for method in ['design','breadth','greedy']:
                chosen=row[f'{method}_phages'].split('|')
                assert len(chosen)==len(set(chosen))==3
                observed=[lookup.loc[(str(row['host']),p)] for p in chosen]
                assert int(any(v==1 for v in observed))==row[f'{method}_lower']==row[f'{method}_upper']
        summary_row=next(r for r in report['summaries'] if r['bank']==bank and r['grouping']==grouping)
        for method in ['design','breadth','greedy']:
            assert np.isclose(frame[f'{method}_lower'].mean(),summary_row[f'{method}_coverage'])
    print('PASS: 60 canonical fits, identical 398-row test cohort, matched training counts, recomputed discrimination, and every additional-bank selected outcome')


if __name__=='__main__': main()
