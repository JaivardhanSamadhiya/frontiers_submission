"""Post hoc stronger baseline on the preserved nested Picard evaluation folds."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/results_v2/post_rejection_audit'
spec=importlib.util.spec_from_file_location('nested',ROOT/'experiments/36_nested_corrected_selector.py')
nested=importlib.util.module_from_spec(spec); spec.loader.exec_module(nested)


def main():
    old=pd.read_csv(OUT/'corrected_nested_selector_hosts.csv')
    up=ROOT/'external/gaborieau_coli_phage_interactions_2023'
    mat=pd.read_csv(up/'data/interactions/interaction_matrix.csv',sep=';',index_col=0)
    meta=pd.read_csv(up/'data/genomics/phages/guelin_collection.csv',sep=';').set_index('phage').loc[mat.columns]
    y=(mat.fillna(0).to_numpy()>0).astype(int); observed=mat.notna().to_numpy()
    names={n:i for i,n in enumerate(mat.index)}
    genus=meta.Genus.fillna('Unknown').to_numpy(); isolation=meta.Phage_host.fillna('Unknown').to_numpy()
    rows=[]; summary=[]
    for threshold,frame in old.groupby('threshold'):
        local=[]
        for fold,test in frame.groupby('fold'):
            train=frame[frame.fold!=fold]
            assert set(train.group).isdisjoint(test.group)
            tr=np.array([names[n] for n in train.bacteria])
            pos=(y[tr]>0)&observed[tr]; covered=np.zeros(len(tr),bool); selected=[]
            prevalence=np.divide(pos.sum(axis=0),observed[tr].sum(axis=0),out=np.zeros(y.shape[1]),where=observed[tr].sum(axis=0)>0)
            for _ in range(3):
                candidates=[j for j in range(y.shape[1]) if j not in selected and not any(genus[j]==genus[a] and isolation[j]==isolation[a] for a in selected)]
                chosen=max(candidates,key=lambda j:(int((pos[:,j]&~covered).sum()),prevalence[j],-j))
                selected.append(chosen); covered|=pos[:,chosen]
            for r in test.to_dict('records'):
                lo,hi,frac,complete=nested.outcome(names[r['bacteria']],selected,y,observed)
                r.update(baseline_lower=lo,baseline_upper=hi,baseline_phages='|'.join(mat.columns[selected]),baseline_complete=complete)
                local.append(r)
        result=pd.DataFrame(local)
        rows.extend(local)
        summary.append({'threshold':threshold,'hosts':len(result),'design_lower':float(result.design_lower.mean()),
            'design_upper':float(result.design_upper.mean()),'greedy_lower':float(result.baseline_lower.mean()),
            'greedy_upper':float(result.baseline_upper.mean()),'paired_bootstrap':nested.paired_bootstrap(result)})
    pd.DataFrame(rows).to_csv(OUT/'picard_greedy_baseline_hosts.csv',index=False)
    (OUT/'picard_greedy_baseline_summary.json').write_text(json.dumps({'status':'post hoc stronger baseline; original nested predictions unchanged',
        'rule':'greedy observed training coverage with same genus AND isolation-host exclusion and three phages',
        'summaries':summary},indent=2))
    print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
