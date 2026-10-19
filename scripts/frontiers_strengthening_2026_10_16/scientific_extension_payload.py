"""Stage a strictly allowlisted scientific extension; never publish or commit.

Paths in the overlay preserve the repository's scientific script/data layout.
The earlier DOI identifies a prior snapshot, not these new analyses.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Run without Python optimization; verification assertions are mandatory')

ROOT=Path(__file__).resolve().parents[2]
STEM=Path('scripts/frontiers_strengthening_2026_10_16')
DATA=Path('data/results_v2/frontiers_strengthening_2026_10_16')
OUT=ROOT/'output/frontiers_strengthening_2026_10_17'
STAGE=OUT/'scientific_extension_staging'
GABORIEAU=Path('external/gaborieau_coli_phage_interactions_2023')
GPIN='99656d8c020fcd1dbbc30a02ee10274a07655562'
GP=DATA/'comparator/official_full_budget_resumed'
PRIOR='7937ee3e705609c2759b37973f9c0d2e97b74f33'
MAX_FILE=40*1024**2
MAX_TOTAL=150*1024**2
CODE=['model_common.py','secondsource_matched.py','model_nested_comparison.py',
      'model_verify.py','model_report.py','verify_extension_analysis.py','doublet_validation.py','test_doublet_validation.py',
      'verify_strengthening.py','verify_comparator_metrics.py','comparator_public_replay.py','comparator_resume_genophi.py','scientific_extension_payload.py','plot_strengthening_doublets.py',
      'model_mlp_budget_sensitivity.py','model_mlp_budget_before_verifier_repair.py','verify_mlp_budget_sensitivity.py','MLP_BUDGET_AMENDMENT.md',
      'DOUBLET_PROTOCOL.md','comparator_external_cohort_audit.md']
DISALLOWED=('manuscript','response_to','reviewer','ai_disclosure','ai_assistance','private','work_plan',
            'build_strengthened','keep_analysis_awake','run_overnight')
TEXT_SUFFIX={'.py','.md','.json','.csv','.txt','.yml','.yaml','.toml','.cff'}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def json_bytes(obj):return (json.dumps(obj,indent=2,ensure_ascii=False)+'\n').encode('utf-8')
def relative(p):
    p=Path(str(p).replace('\\','/'))
    if p.is_absolute() or '..' in p.parts:raise ValueError(f'Nonrelative payload path: {p}')
    return p

def source_pin():
    result=subprocess.run(['git','-c',f'safe.directory={(ROOT/GABORIEAU).as_posix()}',
        '-C',str(ROOT/GABORIEAU),'rev-parse','HEAD'],check=True,capture_output=True,text=True)
    if result.stdout.strip()!=GPIN:raise RuntimeError('Unexpected Gaborieau source commit')
    rows=[]
    for name in ['dev/cocktails/data/cml_cocktails.csv','dev/cocktails/data/cocktails_composition.csv',
                 'data/interactions/interaction_matrix.csv','README.md']:
        blob=subprocess.run(['git','-c',f'safe.directory={(ROOT/GABORIEAU).as_posix()}',
            '-C',str(ROOT/GABORIEAU),'show',f'{GPIN}:{name}'],check=True,capture_output=True).stdout
        disk=(ROOT/GABORIEAU/name).read_bytes()
        # Git may normalize checkout newlines; retain both checksums explicitly.
        if blob.replace(b'\r\n',b'\n')!=disk.replace(b'\r\n',b'\n'):
            raise RuntimeError(f'Source differs from pinned commit: {name}')
        rows.append({'path':(GABORIEAU/name).as_posix(),'sha256':hashlib.sha256(disk).hexdigest(),
                     'git_blob_sha256':hashlib.sha256(blob).hexdigest(),
                     'url':f'https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{GPIN}/{name}'})
    return rows

def comparator_gate():
    folder=ROOT/GP;ap=folder/'audit.json';cp=ROOT/DATA/'comparator/completion.json'
    pending={'state':'pending','full_budget_comparison_complete':False,
             'scope':'Adapted GenoPHI on five Pseudomonas held-out-host folds; no completed comparative claim until audit passes.',
             'reason':'All five outer predictions and a passing independent audit are required.'}
    if not ap.exists() or not cp.exists():
        pending['prediction_locks_present']=[i for i in range(5) if (folder/f'0.001_fold{i}/prediction_lock.json').exists()]
        return pending,[]
    audit=read_json(ap);complete=read_json(cp)
    if complete.get('state')!='all_five_folds_complete_and_audited' or not audit.get('passed'):
        return pending,[]
    assert audit['folds']==5 and audit['selected_models_reloaded']==250 and audit['pairs']==437
    assert len(audit['model_sha256'])==250
    assert len(audit['training'])==5
    assert all(r['selection_runs']==25 and r['selected_models_reloaded']==50 and r['eligible_cutoffs']>0 for r in audit['training'])
    assert audit['all_model_searches_checked']==sum(50*r['eligible_cutoffs'] for r in audit['training'])
    assert audit['auditor_sha256']==sha(ROOT/'experiments/49_audit_official_genophi_outer.py')
    for name,digest in audit['results_sha256'].items():assert sha(folder/relative(name))==digest
    for name,digest in audit['model_sha256'].items():assert sha(ROOT/relative(name))==digest
    paths=[GP/n for n in ['audit.json','matched_metrics.csv','verified_predictions.csv','screening_hosts.csv',
        'model_search_audit.json','GENOPHI_OFFICIAL_FULL_BUDGET_RESULTS.md','run_manifest.json',
        'execution_amendment.json','negative_mcc_repair_manifest.json','driver_before_prediction_api_fix.py']]
    for fold in range(5):
        work=folder/f'0.001_fold{fold}';lock=read_json(work/'prediction_lock.json')
        assert sha(work/'outer_predictions.csv')==lock['prediction_sha256']
        meta=read_json(work/'inputs.json')
        development=relative(meta['development_input'])
        assert sha(ROOT/development)==meta['development_input_sha256']
        paths.append(development)
        paths.append(GP/f'0.001_fold{fold}/selection/features_occurrence.csv')
        paths.append(GP/f'0.001_fold{fold}/development_predictions.csv')
        for pattern in ['selection/run_*/feature_importances.csv','tables/*.csv',
                        'modeling/cutoff_*/run_*/model_performance.csv',
                        'development_summary/model_performance/*.csv']:
            paths.extend(p.relative_to(ROOT) for p in work.glob(pattern) if p.stat().st_size>2)
        for name in ['inputs.json','outer_query_features.csv','development_labels.csv','phage_development_features.csv',
                     'strain_development_features.csv','phage_feature_definitions.csv','strain_feature_definitions.csv',
                     'outer_predictions.csv','prediction_lock.json']:
            paths.append(GP/f'0.001_fold{fold}'/name)
    return {'state':'completed_and_audited','full_budget_comparison_complete':True,'audit_sha256':sha(ap),
        'folds':5,'pairs':437,'hosts':23,'phages':19,'selected_models_verified_locally':250,
        'distribution_scope':'Compact predictions, metrics, input tables, source/selection provenance and model hashes. Candidate trees and 250 selected comparator fitted models are not redistributed; included GBM/MLP outer pipelines are enumerated in the extension manifest.',
        'scope':audit['scope']},paths

def scan(path,content):
    lower=path.as_posix().lower()
    if any(s in lower for s in DISALLOWED) or path.suffix.lower() in {'.docx','.doc','.log','.ps1'}:
        raise ValueError(f'Disallowed file in scientific payload: {path}')
    if len(content)>MAX_FILE:raise ValueError(f'Oversized file: {path}')
    if path.suffix.lower() in TEXT_SUFFIX:
        text=content.decode('utf-8-sig')
        patterns=[r'(?i)[a-z]:[\\/](?:users|documents and settings)[\\/][^\s"<>]+',
                  r'gh[pousr]_[A-Za-z0-9]{30,}',r'github_pat_[A-Za-z0-9_]{30,}',
                  r'sk-[A-Za-z0-9]{40,}',r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----']
        for pattern in patterns:
            if re.search(pattern,text):raise ValueError(f'Private-path/credential-shaped content in {path}')

def check_staged_science():
    """Replay the actual staged payload without changing its result artifacts."""
    code="""
import sys,json,io,contextlib,unittest
sys.path.insert(0,sys.argv[1])
import model_verify as model
import verify_strengthening as independent
with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
    record=model.main(write=False,include_toy_fit=False)
assert record.get('passed') is True
suite=unittest.defaultTestLoader.discover(sys.argv[1],pattern='test_doublet_validation.py')
test_output=io.StringIO()
result=unittest.TextTestRunner(stream=test_output).run(suite)
assert result.wasSuccessful(),test_output.getvalue()
checked={'models':record,'doublet_unit_tests':result.testsRun,
 'independent_doublets':independent.audit_doublets(),
 'independent_secondsource':independent.audit_models(independent.DATA/'models/secondsource'),
 'independent_nested':independent.audit_models(independent.DATA/'models/nested_comparison',nested=True),
 'passed':True,'source':'Checks executed from staged code and staged data; no refitting.'}
if (independent.DATA/'models/mlp_budget_sensitivity').exists():
    import verify_mlp_budget_sensitivity
    checked['mlp_budget_sensitivity']=verify_mlp_budget_sensitivity.verify(write=False)
if (independent.DATA/'comparator/official_full_budget_resumed/audit.json').exists():
    import verify_comparator_metrics
    checked['comparator_metrics']=verify_comparator_metrics.verify()
import pytest
repair_output=io.StringIO()
with contextlib.redirect_stdout(repair_output),contextlib.redirect_stderr(repair_output):
    repair_result=pytest.main([str(independent.ROOT/'tests/test_genophi_retention.py'),'-q','-p','no:cacheprovider'])
assert repair_result==0,repair_output.getvalue()
checked['negative_mcc_repair_tests_passed']=True
print(json.dumps(checked))
"""
    process=subprocess.run([sys.executable,'-B','-c',code,str(STAGE/STEM)],
        check=False,capture_output=True,text=True,timeout=900)
    if process.returncode:raise RuntimeError('Staged scientific checks failed: '+process.stderr)
    return json.loads(process.stdout)

def readme(status):
    d=read_json(ROOT/DATA/'doublets/summary.json')
    evidence=read_json(ROOT/DATA/'models/model_evidence_summary.json')
    ss={r['regime']:r for r in evidence['secondsource'] if r['method']=='full'}
    shared,dual=ss['shared_entities_allowed'],ss['both_clusters_excluded']
    contrast=next(r for r in d['primary_decision'] if r['resampling']=='host_cluster' and r['metric']=='gain_vs_uniform')
    return f'''# PrecisionPhage scientific extension: October 2026

This local staging tree contains completed retrospective scientific additions and their reproducibility inputs. It preserves repository-relative `scripts/`, `data/`, `src/` and `external/` paths. It is an overlay for the existing PrecisionPhage scientific repository, not a manuscript submission bundle. No publication or Git write is performed by the staging command.

The earlier release `v2026.10.15.2`, commit `{PRIOR}`, is archived at https://doi.org/10.5281/zenodo.23426457. **That DOI does not cover these new files.** A new release and version DOI are required when this extension is published. Do not cite an unpublished staging tree as an archived release.

## Completed findings

- Measured equal-size menus: 300 doublets across 100 hosts, with 224 active assays and 26 hosts whose three doublets differ in binary activity. The repaired score selected an active option for 78/100 hosts, compared with expected uniform activity {d['uniform_choice_activity']:.6f}. The gain was {contrast['estimate']:.6f} (paired host-cluster 95% interval {contrast['ci_low']:.6f} to {contrast['ci_high']:.6f}); it does not establish an improvement. The parent triplets were externally selected; this is not arbitrary cocktail design or a new independent cohort.
- Fixed-composition host discrimination: only three of 18 compositions contain both classes. Repaired conditional pair-weighted AUROC is 0.720930; the host-cluster interval is 0.522914–0.875000, but the crossed host/composition interval is 0.250000–0.940000. Shared phages remain part of a fixed bank, not independent bank replication.
- StaphStudy second-source analysis: 1,053 input pairs, canonical 8 phage/26 host clusters, four eligible blocks and 200 evaluated test pairs. Full-model mean pooled AUROC changes from {shared['auroc_mean']:.6f} with shared entities to {dual['auroc_mean']:.6f} with dual exclusion; macro-fold AUROC changes from {shared['macro_fold_auroc_mean']:.6f} to {dual['macro_fold_auroc_mean']:.6f}. Intermediate regimes are not monotonically ordered. Phage-only pooled and within-fold trends disagree, and the training-prior control has pooled AUROC 0.322495 despite within-fold AUROC 0.500. This is a retrospective within-source refit, not new frozen cross-source validation.
- Fresh nested comparison with identical 448-feature inputs: GBM/MLP host-exclusion AUROC 0.935105/0.891437 and dual-exclusion AUROC 0.763409/0.723116. Three candidates per family and identical grouped inner splits were used. All 100 MLP fits reached the declared 150-epoch budget; optimal convergence or exhaustive search is not claimed. This does not reconstruct the historical model-selection rationale.

## Comparator state

GenoPHI state at staging: **{status['state']}**. `COMPARATOR_PUBLIC_STATUS.json` records the gate. If pending, no full-budget comparator result is included or asserted. A completed state requires five outer folds, successful independent representation/model-reload checks, and matching output hashes. Even after completion, this is an adapted Pseudomonas fixed-bank comparison, not an exact end-to-end published GenoPHI reproduction and not the missing all-100 Gaborieau comparator.

The complete candidate model forest, runtime logs, full genomes and private correspondence are excluded. If comparator results are included, its compact tables, prediction locks, audit and model checksums are supplied; its 250 fitted ensemble models are not redistributed here. The 164 fitted GBM/MLP models for the completed second-source and nested analyses are included.

## Reproduction and verification

Use a separate checkout of the prior scientific release and overlay this tree, preserving paths. The science input dependencies used by the quick checks are included, together with the new code and scientific `src/precisionphage` modules. Python 3.11 and recorded package versions are used. `SCIENTIFIC_EXTENSION_ENVIRONMENT.json` records the packaging/audit environment; run-specific protocols retain the training settings and versions where recorded.

From that checkout:

```text
python scripts/frontiers_strengthening_2026_10_16/scientific_extension_payload.py --verify-stage --stage-root . --allow-additional-files
python -B scripts/frontiers_strengthening_2026_10_16/verify_extension_analysis.py
python -B -m unittest discover -s scripts/frontiers_strengthening_2026_10_16 -p test_doublet_validation.py -v
```

The strict staging verifier runs against its default local `output/frontiers_strengthening_2026_10_17/scientific_extension_staging` folder. For an extracted extension-only tree use `--verify-stage --stage-root .` before overlaying it. In a complete repository use the explicit `--allow-additional-files` option shown above: it verifies every extension hash but does not assess other base-repository files, which require their own release checks. The staging verifier writes only its verification report. The read-only `verify_extension_analysis.py` command recalculates controls/metrics and reloads 174 supplied models without fitting or writing reports; it includes 88 exclusion checks, eight original input hashes, candidate selection and the later MLP sensitivity. The separate synthetic two-tree scaler-guard test was run during development and can be rerun with `model_verify.py --read-only`, which fits only that toy example, not any analyzed model. Direct legacy `model_verify.py`, `verify_strengthening.py` and `verify_mlp_budget_sensitivity.py` command-line calls regenerate reports; use the new wrapper to preserve archive bytes. Prediction CSVs must use `float_precision='round_trip'` when recomputing tied-score rank metrics.

To rerun the measured-menu analysis or bounded nested fits, use `doublet_validation.py` and `model_nested_comparison.py`. To reconstruct StaphStudy descriptors and rerun `secondsource_matched.py`, first restore the genome sources described by the earlier release DATA.md; full FASTAs are deliberately absent. Saved descriptor arrays allow the quick model verification without those genomes. Training runs write only their dated analysis output directories. Regeneration may alter runtime/provenance timestamps; preserve release copies when checking exact artifact hashes.

For the full-budget comparator, read `GENOPHI_COMPARISON_PROTOCOL.md`, `GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md` and the dated repair amendment. Use a separate environment with `requirements-genophi-comparator.txt`. Input acquisition/preparation is supplied through `fetch_additional_bank_inputs.py`, `prepare_genophi_comparison.py`, `prepare_official_genophi.py`, `check_genophi_runtime.py` and `verify_genophi_preparation.py` under `scripts/`; manifests identify exact source/protein checksums. The October 5 reduced-budget amendment is historical context, not the full-budget recipe. The Windows runtime smoke helper is platform-specific.

`comparator_resume_genophi.py`, the historical driver and the full training auditor require the original prepared checkpoint tree, repair-chain snapshots, protein caches and fitted candidate models. In particular, the one-time repair registration validates 82 inherited models; it is not a clean-install setup command. Those large historical artifacts are deliberately not in this compact overlay. Do not present the resume commands as a self-contained from-scratch training reproduction. The independent local audit verifies those retained objects; the compact public replay below has the narrower scope stated explicitly.

When the comparator is marked complete, `python scripts/frontiers_strengthening_2026_10_16/verify_comparator_metrics.py` replays its compact saved metrics and paired host-group intervals, checks input-table hashes and prediction alignment, and starts no model fits. The per-fold merged development matrices, feature-occurrence/importance tables, eligible cutoff tables and every completed eight-candidate search table accompany that completed result. This metric-only replay is separate from the full local 250-model reload and protein-representation audit; the latter still requires the deliberately omitted large caches and fitted comparator models.

If a separately identified fitted-audit companion is supplied with the archived release, extract it into the same core directory and use the comparator environment to run `python scripts/frontiers_strengthening_2026_10_16/comparator_public_replay.py`. Its default requires all five completed folds. It verifies retained search models, selected-ensemble predictions, source hashes, development cutoff selection and cached protein-family representations without fitting. Historical forest paths are resolved to byte-identical resumed copies; path-separator normalization is explicit. This narrower public scientific replay intentionally does not certify private execution-log history or reconstruct raw-proteome annotation/alignment/clustering. Explicit `--folds` arguments certify only those partial folds, never the complete comparison. Do not infer that a companion has been published merely because its replay code is present; check the release file inventory and verification record.

## Source and limits

`SCIENTIFIC_EXTENSION_SOURCE_NOTICES.md` and `SCIENTIFIC_EXTENSION_SOURCE_PROVENANCE.json` identify original data, source commits and retained terms. The Gaborieau assay tables are unchanged pinned source files. Derived results retain their original paths. The scientific cohort-search report found no currently eligible untouched cohort for the repaired Picard fixed-bank scorer; no positive transfer claim is manufactured from existing observations.

No manuscript drafts, reviewer-response text, writing scripts, raw AI logs, private mail, user-profile paths or credentials are intended in this scientific payload. `SCIENTIFIC_EXTENSION_MANIFEST.json` inventories all payload hashes; `SCIENTIFIC_EXTENSION_VERIFICATION.json` records the local stage check.
'''

def stage_payload():
    OUT.mkdir(parents=True,exist_ok=True)
    STAGE.mkdir(exist_ok=True)
    assert STAGE.resolve().is_relative_to(OUT.resolve())
    checks=read_json(ROOT/DATA/'models/model_verification.json')
    independent=read_json(OUT/'INDEPENDENT_ANALYSIS_CHECKS.json')
    assert checks['passed'] and checks['artifact_reload_predictions']==164 and independent['all_three_complete']
    sensitivity=ROOT/DATA/'models/mlp_budget_sensitivity'
    additional_models=0
    if sensitivity.exists():
        if not (sensitivity/'completion.json').exists() or not (sensitivity/'verification.json').exists():
            raise RuntimeError('MLP training-budget sensitivity is still pending; do not stage partial fits.')
        sensitivity_checks=read_json(sensitivity/'verification.json')
        assert sensitivity_checks['passed'] and sensitivity_checks['outer_models_reloaded']==10
        assert sha(sensitivity/'completion.json')==sensitivity_checks['completion_sha256']
        additional_models=10
    sources=source_pin();status,comparator_files=comparator_gate()
    paths=set(STEM/name for name in CODE)
    for p in (ROOT/DATA/'models').rglob('*'):
        if p.is_file() and p.name!='MANUSCRIPT_INSERTS.md' and p.suffix.lower() in {'.csv','.json','.npz','.joblib','.md'}:paths.add(p.relative_to(ROOT))
    for p in (ROOT/DATA/'doublets').iterdir():
        if p.is_file() and p.name!='MANUSCRIPT_INSERTS.md' and p.suffix.lower() in {'.csv','.json','.npz','.pdf','.png','.tiff','.md'}:paths.add(p.relative_to(ROOT))
    # Every compact dependency used by the saved scientific audits is explicit.
    for record in read_json(ROOT/DATA/'doublets/summary.json')['inputs'].values():
        p=relative(record['path']);assert sha(ROOT/p)==record['sha256'];paths.add(p)
    for analysis in ['secondsource','nested_comparison']:
        for name,digest in read_json(ROOT/DATA/'models'/analysis/'retrospective_protocol.json')['inputs'].items():
            p=relative(name);assert sha(ROOT/p)==digest;paths.add(p)
    paths.update([Path('data/results_v2/external_staph_predictions.csv'),Path('data/results_v2/external_validation_models/full.joblib'),
        Path('experiments/35_matched_cohort_leakage_audit.py'),Path('experiments/47_official_genophi_outer_comparison.py'),
        Path('experiments/49_audit_official_genophi_outer.py'),Path('GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md'),
        Path('LICENSE'),GABORIEAU/'README.md',DATA/'comparator/source_data_license.json'])
    paths.update(Path(name) for name in ['requirements-genophi-comparator.txt',
        'scripts/prepare_genophi_comparison.py','scripts/prepare_official_genophi.py',
        'scripts/verify_genophi_preparation.py','scripts/check_genophi_runtime.py',
        'scripts/fetch_additional_bank_inputs.py','GENOPHI_COMPARISON_PROTOCOL.md',
        'GENOPHI_TRAINING_AMENDMENT_2026_10_05.md',
        'scripts/finish_official_genophi_comparison.py','scripts/register_negative_mcc_repair.py',
        'tests/test_genophi_retention.py','GENOPHI_NEGATIVE_MCC_AMENDMENT_2026_10_06.md'])
    paths.update(p.relative_to(ROOT) for p in (ROOT/'src/precisionphage').rglob('*.py') if '__pycache__' not in p.parts)
    for name in ['official_source_manifest.json','matched_outer_folds.csv','preparation_checks.json','upstream_source.json',
                 'input_manifest.json','runtime_check.json']:
        paths.add(Path('data/results_v2/genophi_comparison')/name)
    official=read_json(ROOT/'data/results_v2/genophi_comparison/official_source_manifest.json')
    for entry in official['files']:
        p=Path('data/external_validation/genophi_comparison/official_source')/relative(entry['path'])
        assert sha(ROOT/p)==entry['sha256'];paths.add(p)
    paths.update(comparator_files)
    rows=[]
    for p in sorted(paths):
        source=ROOT/p
        if not source.is_file():raise FileNotFoundError(source)
        content=source.read_bytes();scan(p,content)
        dest=STAGE/p;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
        assert sha(source)==sha(dest)
        rows.append({'path':p.as_posix(),'bytes':len(content),'sha256':sha(dest),'provenance':'byte-identical local scientific source'})
    scientific_checks=check_staged_science()
    notices='''# Scientific extension source notices

PrecisionPhage code retains the MIT license. Third-party data and software are not relicensed.

Gaborieau et al. (2024), Nature Microbiology 9:2847–2861, https://doi.org/10.1038/s41564-024-01832-5. Original assay/composition/Picard tables and public README are pinned to commit 99656d8c020fcd1dbbc30a02ee10274a07655562 at https://github.com/mdmparis/coli_phage_interactions_2023. The archived code/data record https://zenodo.org/records/13831957 identifies Creative Commons Attribution 4.0, https://creativecommons.org/licenses/by/4.0/. Attribution and original terms continue to apply. No publication PDF or private correspondence is redistributed.

VHIP source: Bastien et al. (2024), https://doi.org/10.1371/journal.pcbi.1011649; source experiment StaphStudy: Goller et al. (2021), https://doi.org/10.1038/s41467-021-27037-6. The previous release provides full accession and source provenance. The unchanged compact input/cache artifacts accompany these new retrospective fits.

GenoPHI: Noonan et al. (2026), https://doi.org/10.1038/s41564-026-02482-5. Official scientific source files are pinned to f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f and included byte-identically with their LICENSE at data/external_validation/genophi_comparison/official_source/LICENSE. The source manifest provides every raw URL/hash. The processed public data record https://doi.org/10.5281/zenodo.21696332 identifies Creative Commons Attribution 4.0; its original attribution and terms apply. The dated source_data_license.json records the original API metadata and archive checksum. Execution adaptations, including negative-MCC retention, are documented separately; this is not an exact end-to-end published reproduction. Large genomes, protein caches and candidate forests are omitted from this core overlay.

Only scientific code, original compact data with attribution, derived results, provenance and verification are staged. The prior DOI 10.5281/zenodo.23426457 does not identify this extension.
'''
    readme_text=readme(status)
    if additional_models:
        readme_text=readme_text.replace('The 164 fitted GBM/MLP models for the completed second-source and nested analyses are included.',
            'The 164 fitted GBM/MLP models for the completed second-source and nested analyses and ten additional MLP budget-sensitivity outer models are included (174 total).')
        sensitivity_metrics=__import__('pandas').read_csv(sensitivity/'pooled_metrics.csv',float_precision='round_trip')
        extended=sensitivity_metrics.loc[sensitivity_metrics.method.eq('MLP_extended')].set_index('regime')
        sensitivity_text=(f"## MLP optimization-budget sensitivity\n\nA prewritten retrospective amendment retained the same grouped splits and three neural candidates, "
            f"allowing at most 1,000 epochs with training-loss stopping (no random validation). MLP AUROC was "
            f"{extended.loc['host_cluster','auroc']:.6f} for host exclusion and {extended.loc['dual_cluster','auroc']:.6f} for dual exclusion. "
            f"The 10 additional outer models and all traces are included. Observed epochs ranged from "
            f"{sensitivity_checks['epochs_min']} to {sensitivity_checks['epochs_max']}; {sensitivity_checks['fits_reaching_1000']} of 100 fits reached the new cap. "
            "The initial 150-epoch results remain unchanged. A documented CSV-verification tolerance repair changed no fit or score. "
            "This is an optimization sensitivity on already accessed outcomes, not exhaustive tuning or fresh external validation.\n\n")
        readme_text=readme_text.replace('## Comparator state',sensitivity_text+'## Comparator state')
    generated={
        Path('SCIENTIFIC_EXTENSION_README.md'):readme_text.encode('utf-8'),
        Path('SCIENTIFIC_EXTENSION_SOURCE_NOTICES.md'):notices.encode('utf-8'),
        Path('SCIENTIFIC_EXTENSION_SOURCE_PROVENANCE.json'):json_bytes({'gaborieau_commit':GPIN,'source_files':sources,'genophi_commit':official['commit']}),
        Path('COMPARATOR_PUBLIC_STATUS.json'):json_bytes(status),
        Path('SCIENTIFIC_EXTENSION_SCIENTIFIC_CHECKS.json'):json_bytes(scientific_checks),
        Path('SCIENTIFIC_EXTENSION_ENVIRONMENT.json'):json_bytes({'scope':'Packaging/verification runtime; run protocols retain their own settings.',
            'packages':{name:importlib.metadata.version(name) for name in ['numpy','pandas','scikit-learn','xgboost','scipy','joblib','matplotlib','sourmash']}}),
        Path('reproducibility/strengthening_2026_10_17/INDEPENDENT_ANALYSIS_CHECKS.json'):json_bytes(independent)}
    for p,content in generated.items():
        scan(p,content);dest=STAGE/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(content)
        rows.append({'path':p.as_posix(),'bytes':len(content),'sha256':sha(dest),'provenance':'generated scientific packaging metadata'})
    if sum(r['bytes'] for r in rows)>MAX_TOTAL:raise ValueError('Payload exceeds compact staging budget')
    expected={r['path'] for r in rows}|{'SCIENTIFIC_EXTENSION_MANIFEST.json','SCIENTIFIC_EXTENSION_VERIFICATION.json'}
    extras=[p.relative_to(STAGE).as_posix() for p in STAGE.rglob('*') if p.is_file() and p.relative_to(STAGE).as_posix() not in expected]
    if extras:raise RuntimeError(f'Unexpected pre-existing staging files; no deletion performed: {extras[:20]}')
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'prior_release':'v2026.10.15.2','prior_commit':PRIOR,
        'prior_doi':'10.5281/zenodo.23426457','new_release_or_doi_assigned':False,'comparator':status,
        'copied_scientific_outer_models':164+additional_models,'files':sorted(rows,key=lambda r:r['path']),
        'total_bytes':sum(r['bytes'] for r in rows),'excluded':['manuscript and reviewer material','writing scripts','runtime/AI/mail logs',
            'full genomes/protein caches','GenoPHI candidate trees and comparator model binaries','user-profile paths/credentials'],
        'layout':'Overlay preserves scientific source and analysis data paths; quick verification dependencies included.'}
    (STAGE/'SCIENTIFIC_EXTENSION_MANIFEST.json').write_bytes(json_bytes(manifest))
    return verify_stage(STAGE)

def verify_stage(stage,allow_additional_files=False):
    stage=Path(stage).resolve();manifest=read_json(stage/'SCIENTIFIC_EXTENSION_MANIFEST.json')
    for row in manifest['files']:
        p=relative(row['path']);f=stage/p
        assert f.is_file() and f.stat().st_size==row['bytes'] and sha(f)==row['sha256'],p
        scan(p,f.read_bytes())
    fitted=list((stage/DATA/'models').rglob('*.joblib'));assert len(fitted)==manifest['copied_scientific_outer_models']
    expected={row['path'] for row in manifest['files']}|{'SCIENTIFIC_EXTENSION_MANIFEST.json','SCIENTIFIC_EXTENSION_VERIFICATION.json'}
    actual={p.relative_to(stage).as_posix() for p in stage.rglob('*') if p.is_file()}
    extra=actual-expected
    if not allow_additional_files:assert not extra,sorted(extra)
    report={'passed':True,'verified_utc':datetime.now(timezone.utc).isoformat(),'files':len(manifest['files']),
        'bytes':manifest['total_bytes'],'outer_models':len(fitted),'manifest_sha256':sha(stage/'SCIENTIFIC_EXTENSION_MANIFEST.json'),
        'comparator_complete':manifest['comparator']['full_budget_comparison_complete'],
        'public_write_performed':False,'private_or_unexpected_files_found':False,
        'additional_files_not_assessed':len(extra),
        'scope':'Extension allowlist only; additional repository files are not assessed' if allow_additional_files else 'Strict extension-only tree'}
    (stage/'SCIENTIFIC_EXTENSION_VERIFICATION.json').write_bytes(json_bytes(report))
    print(json.dumps(report,indent=2));return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify-stage',action='store_true')
    parser.add_argument('--stage-root',type=Path,default=STAGE)
    parser.add_argument('--allow-additional-files',action='store_true');args=parser.parse_args()
    if args.verify_stage:verify_stage(args.stage_root,args.allow_additional_files)
    else:
        if args.stage_root!=STAGE or args.allow_additional_files:raise ValueError('Staging destination is fixed; additional-file permission is verification-only')
        stage_payload()
