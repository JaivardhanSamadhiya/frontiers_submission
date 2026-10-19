"""Independently verify proteome provenance and preserve matched outer folds."""
from pathlib import Path
import json
import zipfile
import pandas as pd
from prepare_genophi_comparison import ROOT, OUT, MATRIX, EXPECTED, sha, protein_records, immutable_write


def main():
    manifest=json.loads((OUT/'input_manifest.json').read_text())
    assert manifest['complete'] and not manifest['failures']
    assert sha(MATRIX.read_bytes()) == EXPECTED
    with zipfile.ZipFile(MATRIX) as z:
        data=pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'))
    data=data[data.dataset=='Pseudomonas']
    assert len(data)==437 and data.interaction.sum()==158
    expected={'host':set(data.strain),'phage':set(data.phage)}
    for kind,names in expected.items():
        entries=[e for e in manifest['proteomes'] if e['kind']==kind]
        assert len(entries)==len(names) and {e['entity'] for e in entries}==names
    for entry in manifest['proteomes']:
        raw=(ROOT/entry['path']).read_bytes()
        assert sha(raw)==entry['sha256']
        records=protein_records(raw)
        assert len(records)==entry['proteins']
        assert sum(len(r.seq) for r in records)==entry['amino_acids']
    source=ROOT/'data/results_v2/additional_banks/host_selections.csv'
    predictions=ROOT/'data/results_v2/additional_banks/pair_predictions.csv'
    h=pd.read_csv(source,dtype={'grouping':str})
    p=pd.read_csv(predictions,dtype={'grouping':str})
    folds=h[(h.bank=='Pseudomonas') & h.grouping.isin(['0.001','identity'])][
        ['grouping','host','group','fold']].sort_values(['grouping','host'])
    assert len(folds)==46 and not folds.duplicated(['grouping','host']).any()
    baseline=p[(p.bank=='Pseudomonas') & p.grouping.isin(['0.001','identity'])]
    for grouping,sub in folds.groupby('grouping'):
        assert set(sub.host)==expected['host'] and sub.fold.nunique()==5
        assert sub.groupby('group').fold.nunique().max()==1
        pairs=baseline[baseline.grouping==grouping]
        merged=pairs.merge(data,left_on=['host','phage'],right_on=['strain','phage'],
                           how='outer',validate='one_to_one',indicator=True)
        assert len(merged)==437 and (merged['_merge']=='both').all()
        assert (merged.label==merged.interaction).all()
        joined=pairs.merge(sub,on=['grouping','host'],validate='many_to_one',suffixes=('_old','_locked'))
        assert (joined.fold_old==joined.fold_locked).all()
    foldpath=OUT/'matched_outer_folds.csv'
    immutable_write(foldpath,folds.to_csv(index=False,lineterminator='\n').encode())
    runtime=json.loads((OUT/'runtime_check.json').read_text())
    assert runtime['complete'] and all(r['returncode']==0 for r in runtime['runs'])
    report={'checks_passed':True,'hosts':23,'phages':19,'pairs':437,'positive_pairs':158,
        'proteins':sum(e['proteins'] for e in manifest['proteomes']),
        'matched_grouping_regimes':['0.001','identity'],
        'groups_at_0.001':int(folds[folds.grouping=='0.001'].group.nunique()),
        'fold_assignment_sha256':sha(foldpath.read_bytes()),
        'original_host_results_sha256':sha(source.read_bytes()),
        'original_pair_predictions_sha256':sha(predictions.read_bytes()),
        'input_manifest_sha256':sha((OUT/'input_manifest.json').read_bytes()),
        'stage':'Inputs and matched folds verified; predictive comparison not run.'}
    (OUT/'preparation_checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
