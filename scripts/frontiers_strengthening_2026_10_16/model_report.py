"""Summarize completed, verified model analyses without editing manuscripts."""
import json
import pandas as pd
from model_common import OUT,save_json

def main():
    checks=json.loads((OUT/'model_verification.json').read_text());assert checks['passed']
    s=OUT/'secondsource';n=OUT/'nested_comparison'
    sp=pd.read_csv(s/'pooled_metrics.csv');np=pd.read_csv(n/'pooled_metrics.csv')
    sf=pd.read_csv(s/'fold_metrics.csv');nf=pd.read_csv(n/'outer_metrics.csv')
    sc=json.loads((s/'completion.json').read_text());nc=json.loads((n/'completion.json').read_text())
    audit=json.loads((s/'input_audit.json').read_text());blocks=json.loads((s/'all_block_eligibility.json').read_text())
    summary=sp.groupby(['method','regime']).agg(auroc_mean=('auroc','mean'),auroc_min=('auroc','min'),auroc_max=('auroc','max'),
       AP_mean=('average_precision','mean'),macro_fold_auroc_mean=('macro_fold_auroc','mean'),within_fold_weighted_auroc_mean=('within_fold_pair_weighted_auroc','mean')).reset_index()
    summary.to_csv(s/'seed_summary.csv',index=False)
    text=['# Completed retrospective model strengthening','',
       'These are new retrospective analyses of previously accessed sources. They do not restore the original tuning history, create a new untouched test cohort, or change preserved historical predictions.','',
       '## StaphStudy second-source matched evaluation','',
       f'The input contains {audit["rows"]} pairs ({audit["positives"]} positive), all 39 phages × 27 host reference identifiers. Canonical 21-mer sketches of size 2048 at distance 0.05 produce {audit["phage_clusters"]} phage and {audit["host_clusters"]} host clusters. Four of five predetermined diagonal blocks pass the minimum three-positive/three-negative rule; they contain 200 test pairs (64 positive; prevalence 0.320). Every regime uses those exact test rows and fold-specific matched training class counts.','',
       'Assigned-block eligibility: '+json.dumps(blocks), '',
       'Three sampling seeds are sensitivity runs, not independent datasets. Values below average their pooled and within-fold metrics.','',
       summary.round(6).to_string(index=False),'',
       'Interpretation: the full model worsens from shared to dual exclusion, but the full four-step ordering is not monotonic and phage-only pooled AUROC does not consistently fall. The training-prior control is 0.322 pooled and exactly 0.500 within every fold, demonstrating a cross-fold prevalence artifact. Phage-only pooled and macro-fold trends disagree. Therefore, report pooled, macro-fold and positive-negative-pair-weighted within-fold metrics together; this is evidence of source-dependent evaluation behavior, not a universal hierarchy. Pair-only features outperform the full model under dual exclusion descriptively. Only four eligible blocks limit inference.','',
       f'Runtime: {sc["elapsed_seconds"]:.1f} seconds, {sc["fits"]} CPU fits; 144 fitted tree artifacts retained. Input reconstruction reproduces the historical frozen StaphStudy predictions to maximum absolute error {checks["staph_assembly_matches_frozen_prediction_max_error"]:.3g}.','',
       '## Fresh compatible-input nested GBM versus MLP','',
       'Both families use the same 448 columns, identical outer/inner rows, train-only scaling and three predeclared candidates. Host evaluation uses five group-size-balanced folds; dual evaluation uses the five original dual-exclusion blocks. Inner folds follow the same exclusion axes. Candidate selection uses only mean inner AUROC; there is no random-row fallback, random-row early stopping, calibration or outer-test-based selection. MLP uses 150 fixed epochs; GBM uses 200 or 400 trees. Equal candidate and split counts are not equal wall-clock, capacity, or exhaustive tuning budgets.','',
       np.round(6).to_string(index=False),'',
       'Per-fold results:','',nf[['regime','fold','family','chosen_candidate','n','positive','auroc','average_precision']].round(6).to_string(index=False),'',
       f'Runtime: {nc["elapsed_seconds"]:.1f} seconds; {nc["inner_fits"]} inner fits plus {nc["outer_refits"]} outer refits. {nc["convergence_warning_count"]} recorded warnings reflect fixed-epoch MLP runs reaching their budget; convergence is not claimed. All 20 outer fitted pipelines are retained, with inner scores, exact configurations and all split indices. Inner fitted objects themselves are not retained.','',
       'This answers a new, bounded model-family comparison. It does not prove that the historical GBM was selected by a documented comparison, nor does it establish one family as universally superior. Host and dual outer cohorts differ, so compare families within each regime.','',
       '## Verification','',json.dumps(checks,indent=2),'',
       'No GPU or paid compute was used. No manuscript, portal declaration, public release, or historical result was overwritten by these scripts.']
    (OUT/'MODEL_STRENGTHENING_REPORT.md').write_text('\n'.join(text),encoding='utf-8')
    save_json(OUT/'model_evidence_summary.json',{'secondsource':summary.to_dict('records'),'nested_comparison':np.to_dict('records'),
        'verified':True,'secondsource_inputs':audit,'limits':['retrospective','no universal hierarchy','four second-source blocks','equal candidate counts not compute budgets','old tuning unrecoverable']})

if __name__=='__main__':main()
