"""Independently verify completed frozen/screening analyses, without fitting."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
from math import comb
import joblib
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/results_v2/basel_ecor_validation'
SCREEN = ROOT/'data/results_v2/screening_domain_biology'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def read(path): return json.loads(Path(path).read_text())


def main():
    lock = read(OUT/'prediction_lock.json'); result = read(OUT/'evaluation.json')
    uncertainty = read(OUT/'uncertainty_definition_lock.json')
    for manifest in [lock, read(SCREEN/'manifest.json')]:
        for name, expected in manifest['sha256'].items():
            assert sha(ROOT/name) == expected, name
    assert sha(OUT/'prediction_lock.json') == result['prediction_lock_sha256']
    assert sha(ROOT/'BASEL_ECOR_UNCERTAINTY_AMENDMENT_2026_10_05.md') == uncertainty['protocol_sha256']
    assert sha(OUT/'secondary_host_clusters.csv') == uncertainty['assignments_sha256']
    assert sha(ROOT/'experiments/48_basel_uncertainty_sensitivity.py') == uncertainty['script_sha256']
    assert uncertainty['recorded_utc'] < lock['locked_utc'] < result['evaluated_utc']
    secondary = read(OUT/'secondary_uncertainty.json')
    assert sha(OUT/'evaluation.json') == secondary['evaluation_sha256']
    assert sha(OUT/'uncertainty_definition_lock.json') == secondary['definition_lock_sha256']
    identifier = read(OUT/'identifier_manifest.json')
    raw_path = ROOT/identifier['sealed_outcomes']
    assert sha(raw_path) == identifier['sha256']
    genome_manifest = read(OUT/'genome_download_manifest.json')
    for r in genome_manifest['genomes']:
        fasta = ROOT/'data/external_validation/basel_ecor_2026'/f'{r["kind"]}_genomes'/f'{r["identifier"]}.fasta'
        assert sha(fasta) == r['sha256'], str(fasta)
    assert len(genome_manifest['genomes']) == 77
    pred = pd.read_csv(OUT/'frozen_predictions.csv')
    labels = pd.read_csv(raw_path).rename(columns={'strain': 'host', 'interaction': 'label'})
    merged = pred.merge(labels, on=['host', 'phage'], how='outer', validate='one_to_one')
    assert len(merged) == len(pred) == len(labels) == 1240
    assert not merged[['score', 'label']].isna().any().any()
    scored = pd.read_csv(OUT/'evaluated_predictions.csv')
    alignment = merged.merge(scored, on=['host', 'phage'], validate='one_to_one', suffixes=('', '_saved'))
    np.testing.assert_array_equal(alignment.label, alignment.label_saved)
    np.testing.assert_array_equal(alignment.score, alignment.score_saved)
    frozen = joblib.load(ROOT/'data/results_v2/external_validation_models/sequence_only.joblib')
    X = np.load(OUT/'frozen_input_matrix.npz', allow_pickle=False)['X']
    assert X.shape == (1240, 444) and np.isfinite(X).all()
    restored = frozen['model'].predict_proba(frozen['scaler'].transform(X))[:, 1]
    max_difference = float(np.max(np.abs(restored-pred.score)))
    np.testing.assert_allclose(restored, pred.score, rtol=0, atol=1e-7)
    y, s = merged.label.to_numpy(), merged.score.to_numpy()
    npos, nneg = int(y.sum()), int((y == 0).sum())
    auc = (rankdata(s)[y == 1].sum()-npos*(npos+1)/2)/(npos*nneg)
    brier = float(np.mean((s-y)**2))
    bins = np.minimum((s*10).astype(int), 9)
    ece = sum(abs(float(np.mean(s[bins == b])-np.mean(y[bins == b])))*np.mean(bins == b)
              for b in np.unique(bins))
    np.testing.assert_allclose([auc, average_precision_score(y,s), brier, ece],
                              [result['metrics']['roc_auc'], result['metrics']['pr_auc'],
                               result['metrics']['brier'], result['calibration']['ece']], rtol=0, atol=1e-12)
    assert npos == 302 and nneg == 938 and (s > .5).all()
    assert scored.host_cluster.nunique() == 1 and scored.phage_cluster.nunique() == 31
    assert 'unavailable' in result['ci95']
    training = pd.read_csv(ROOT/'data/interim_v2/interactions_modelable.csv')
    training = training[training.study == 'NCBI_HR']
    ecoli = training[training.host_genus == 'escherichia']
    training_diagnostic = {
        'training_pairs': len(training), 'training_positive': int(training.label.sum()),
        'escherichia_pairs': len(ecoli), 'escherichia_positive': int(ecoli.label.sum()),
        'escherichia_host_identifiers': sorted(ecoli.host.unique().tolist()),
        'new_test_positive_fraction': float(y.mean()), 'new_test_mean_score': float(s.mean()),
        'interpretation': 'Coarse host-reference resolution and prevalence differ; not a proven causal decomposition.'}
    (OUT/'training_resolution_diagnostic.json').write_text(json.dumps(training_diagnostic, indent=2))
    # Independent direct top-k calculation, not the implementation under test.
    original = pd.read_csv(ROOT/'data/results_v2/additional_banks/pair_predictions.csv', dtype={'grouping': str})
    hosts = pd.read_csv(SCREEN/'screening_host_metrics.csv', dtype={'grouping': str})
    checked = 0
    for (bank, grouping, host), z in original.groupby(['bank', 'grouping', 'host']):
        saved = hosts[(hosts.bank == bank) & (hosts.grouping == grouping) & (hosts.host == host)]
        n, positives = len(z), int(z.label.sum())
        for method, col in [('model', 'score'), ('training_prevalence', 'breadth_score'), ('random_expected', None)]:
            ordered = z.sort_values([col, 'phage'], ascending=[False, True]) if col else None
            for r in saved[saved.method == method].itertuples():
                if col:
                    top = int(ordered.iloc[:r.budget].label.sum())
                    expected = [float(top > 0), top/r.budget]
                else:
                    expected = [1-comb(n-positives,r.budget)/comb(n,r.budget) if r.budget <= n-positives else 1, positives/n]
                np.testing.assert_allclose([r.hit,r.precision],expected,rtol=0,atol=1e-12)
                assert r.positive_phages == positives and r.bank_size == n
                if positives == 0: assert np.isnan(r.active_recall_actionable_only) and r.hit == 0
                checked += 1
    assert checked == len(hosts)
    bio = pd.read_csv(SCREEN/'biology_pair_annotations.csv')
    assert len(bio) == 297 and bio.bacteria.nunique() == 99 and 'AN03' not in set(bio.bacteria)
    actual, predicted = bio.productive.to_numpy(), (bio.score_ensemble.to_numpy() >= .5).astype(int)
    counts = {'tp': int(((actual == 1) & (predicted == 1)).sum()),
              'tn': int(((actual == 0) & (predicted == 0)).sum()),
              'fp': int(((actual == 0) & (predicted == 1)).sum()),
              'fn': int(((actual == 1) & (predicted == 0)).sum())}
    assert counts == {'tp':113, 'tn':70, 'fp':60, 'fn':54}
    audit = {'passed': True, 'verified_utc': datetime.now(timezone.utc).isoformat(),
             'frozen_model_reloaded': True, 'maximum_prediction_roundtrip_difference': max_difference,
             'genomes_hash_verified': 77, 'pairs': 1240, 'independent_rank_sum_auroc': float(auc),
             'screening_rows_independently_verified': checked, 'biology_confusion_counts': counts,
             'training_resolution': training_diagnostic,
             'pending': 'Full-budget official GenoPHI comparator is separately audited only after all five folds.',
             'verifier_sha256': sha(Path(__file__)),
             'output_sha256': {str(p.relative_to(ROOT)): sha(p) for folder in [OUT, SCREEN]
                               for p in folder.iterdir() if p.suffix in ['.json', '.csv'] and p.name != 'completed_stage_audit.json'}}
    (OUT/'completed_stage_audit.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps({k:v for k,v in audit.items() if k != 'output_sha256'}, indent=2))


if __name__ == '__main__': main()
