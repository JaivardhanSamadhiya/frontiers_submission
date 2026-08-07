# PrecisionPhage — Independent Scientific & Reproducibility Audit

**Reviewer stance:** Nature Biotechnology / Bioinformatics reviewer, computational virologist,
ML generalization scientist, reproducibility auditor, senior software engineer.
**Scope of this report:** static audit of every source file in `full_pipeline/`.
**Verdict (one line):** *Not ready for peer review.* The code is unusually well-organized and
self-aware about leakage, but it contains at least three publication-fatal methodological
problems and cannot currently be reproduced or executed.

---

## 0. Execution / reproducibility blocker (must read first)

The empirical tasks in the brief (retrain, ablations, generalization benchmarks, SHAP,
figures) **could not be executed**, because the repository as delivered is code-only:

| Requirement | Status |
|---|---|
| `scripts/config.yaml` (every module reads `cfg[...]`) | **MISSING** |
| `data/raw/VirusHostInter.csv` (primary positives) | **MISSING** |
| `data/fastas/hosts/*.fasta`, phage FASTAs | **MISSING** |
| `phagesdb_records_cache.json`, `virushostdb_raw.tsv` | **MISSING** |
| Frozen dataset CSV + `.sha256` | **MISSING** |
| `requirements.txt` / environment lock | **MISSING** |
| `DATA.md` (referenced by `fetch_data.verify_manual_inputs`) | **MISSING** |
| Git repository | **NOT a repo** |
| Python ≥3.10 (code uses `X | None`, f-strings) | only 3.1.1 / 3.4.3 present |
| numpy/pandas/scikit-learn/scipy/matplotlib/torch/PyG | **none installed** |

A reviewer who cloned this would not be able to run anything. This alone is a desk-reject for
a reproducibility-graded venue. Everything below is therefore derived from reading the code;
empirical magnitudes (how much each leak inflates AUC) must be measured once the data + a
working environment are supplied.

---

## CRITICAL TASK 1 — Repository map (data → result)

```
fetch_data.py            verify_manual_inputs(); fetch_all(): VHDB, PhagesDB, INPHARED,
                         NCBI phage FASTAs; ensure_validation_csvs() DERIVES S. aureus
                         "validation" from VirusHostInter.csv   ← see Finding A
        │
data_collection.py       assemble_dataset(): merge VHI + VHDB + INPHARED + NCBI + PhagesDB
                         → POSITIVES ONLY (label=1). _clean_host strips strain/subsp;
                         _clean_phage only lowercases (does NOT strip host tokens ← Finding D)
        │
build_labels.py          build_negatives(): within-genus + cross-genus constructed negatives;
                         q75 imputation of pre-existing numeric cols on negatives ← Finding B
                         freeze_dataset(): SHA-256.  validate_no_leakage(): exact-pair only ← Finding E
        │
data_enrichment.py       compute_all_features(): k-mers, codon, GC, length; PCA(tetra/codon)
                         and SVD(names) FIT ON FULL DATASET ← Finding C; pair features
                         (tetra_corr, cub_dist, gc_match, len_ratio); name + sanitized-name SVD
        │
model.py                 build_feature_matrices(): RAW matrices; per-fold StandardScaler only;
                         _VHI_NUMERIC edge cols [k3dist,k6dist,GCdiff,Homology] ← Finding B
        │
evaluation.py            LOSO over host species; RF / HistGB / XGB; metrics() swallows errors
                         → 0.5/0.0 ← Finding J.  run_experimental_validation() ← Findings A,M
        │
gnn.py                   GAT/SAGE + parallel MLP bypass + learned alpha gate ← Finding G;
                         silent fallback to HistGB on any torch error ← Finding G
        │
logo_validation.py       LOGO over genus.   strain_sim.py: unseen-phage Monte-Carlo.
cocktail_optimizer.py    greedy set-cover on LOSO predictions.
statistical_tests.py     Wilcoxon+BH-FDR, DeLong, McNemar, (cluster) bootstrap ← Finding I
ablation_study.py        feature-group ablation; name-leakage ablation (threshold 0.02) ← Finding D
plotting.py / summary.py figures + methods_skeleton.md + reproducibility_card.md
run_pipeline.py          orchestrates steps 0–9
```

**Train/test split semantics:** LOSO = hold out one cleaned host *species string*; LOGO = hold
out one genus; unseen-phage = remove a random phage set. Negatives are built once globally and
then split, so the *same* feature transforms (PCA/SVD) span folds.

---

## CRITICAL TASK 2 — Leakage audit (assume guilty until proven innocent)

### Finding A (PUBLICATION-FATAL) — "Independent experimental validation" is circular
`fetch_data.ensure_validation_csvs()` (`fetch_data.py:246-289`) builds the S. aureus
"held-out experimental" CSVs **by filtering the very same `VirusHostInter.csv`** that
`data_collection._from_vhi()` (`data_collection.py:116-144`) uses to build training positives.
The docstrings (`data_collection.py:15-17`) and Methods text (`summary.py:141-146`) claim these
pairs are "NEVER used in training" and constitute an "independent experimental validation set."
**They are a subset of the training source.** `run_experimental_validation()` drops only exact
(phage,host) string overlaps (`evaluation.py:294-309`); the remainder are non-overlapping only
due to incidental name-cleaning differences. This is not independent validation and must not be
presented as such.

### Finding B (CRITICAL, latent) — constant-imputed pair features make negatives trivially separable
`build_negatives()` imputes **every pre-existing numeric column** on constructed negatives with
the **75th percentile of the positives** (`build_labels.py:129-144`). `model._VHI_NUMERIC =
[k3dist, k6dist, GCdiff, Homology]` are wired into the edge-feature matrix
(`model.py:30-40, 114-120`). If those VHI-derived columns are present at negative-construction
time (the brief explicitly lists them as model features, so the real/frozen dataset has them),
then **all negatives share one identical constant value** for each while positives vary — a
classifier separates classes from the imputation artifact alone, not biology. In the current
pure-rebuild path `_from_vhi` happens to drop those columns, so the leak is *latent*; it
activates the moment VHI numeric features are ingested or a frozen CSV containing them is used.
This is the single most dangerous landmine in the codebase.
*(Pair features computed in `data_enrichment` — `tetra_corr/cub_dist/gc_match/len_ratio` — are
computed from real sequences for both classes, so they are not subject to this imputation leak.)*

### Finding C (MAJOR) — PCA/SVD fit on the full dataset (preprocessing leakage across folds)
`compute_all_features()` fits PCA on tetranucleotide/codon matrices over **all** phages and
hosts (`data_enrichment.py:436-455`) and SVD on **all** names (`:522-532`) *before* any CV
split. Per-fold `StandardScaler` is correctly fit on train only (`model.py:56-67`), but the
much higher-variance dimensionality reductions are global, so each held-out fold's sequences
shape the representation used to predict them. Unsupervised ⇒ mild, but for a paper whose
headline claim is generalization to unseen entities, reviewers will require these transforms to
be fit inside each training fold.

### Finding D (MAJOR) — name embeddings encode host taxonomy; audit is too lenient
`_clean_phage()` (`data_collection.py:66-72`) does **not** remove host genus/species tokens, so
char-n-gram SVD of phage names (`data_enrichment.py:268-286, 522`) directly encodes the host
(e.g. "staphylococcus phage k"). The sanitization+ablation is good practice, but: (1)
`sanitize_phage_name` only blanks exact host-genus strings present in the host list —
abbreviations, species epithets, and misspellings survive (`data_enrichment.py:289-297`); (2)
the pass criterion "AUC drop < 0.02 ⇒ no leakage" is arbitrary and lenient
(`ablation_study.py:221-227`); (3) the SVD is global (Finding C). **Recommendation:** drop
original-name embeddings from the publication model and report the sanitized/no-name model as
the headline.

### Finding E (MAJOR) — leakage check is exact-pair only; no homology/synonym control
`validate_no_leakage()` compares exact lowercased (phage,host) strings
(`build_labels.py:193-205`). It therefore **cannot** detect: duplicate genomes under different
names, synonymous/alias host names, strain duplicates, or near-identical phages. It "passes"
trivially for LOSO/LOGO and gives false assurance. There is **no sequence-similarity dedup
anywhere** (no MASH/ANI/CD-HIT). LOSO holds out a species string while congeneric near-relatives
remain in training → genus-level homology inflation. LOGO is the right counter-measure and is
implemented; the gap between LOSO and LOGO is the empirical proxy for this leak and must be
reported.

### Finding D/E summary table

| Leakage class | Present? | Mechanism | Mitigation in code | Gap |
|---|---|---|---|---|
| Pair | Guarded | exact-pair check | `validate_no_leakage` | string-only |
| Species | Partial | LOSO holds species string | host cleaning to species | no synonym map |
| Genus | Likely inflates LOSO | congeners remain in train | LOGO provided | report LOSO−LOGO Δ |
| Homology | **Unhandled** | no similarity dedup | none | add MASH/ANI filter |
| Name | **Present** | host tokens in phage names | sanitize + ablation | global SVD, lenient threshold |
| Pair-feature | **Latent-critical** | q75 constant imputation | flag column | drop imputation / recompute |

---

## CRITICAL TASK 3 — Negative-sampling validity

- **False negatives (open-world):** constructed negatives may be unobserved positives. The code
  honestly notes this caveat (`build_labels.py:185-186`) but the metrics treat negatives as
  ground truth. No PU-learning, no confidence weighting.
- **Sampling bias:** cross-genus negatives are almost always genomically dissimilar ⇒ trivially
  easy ⇒ inflate ROC/PR-AUC and create an unrealistic boundary (`build_labels.py:101-122`). The
  `hard_negatives_only` ablation (within-genus only, `ablation_study.py:138-157`) is the correct
  diagnostic and should become a **headline** result, not a side ablation.
- **Class balance:** `n_within * len(pos_hosts)` can explode for promiscuous phages; no global
  prevalence control, so pooled AUC mixes fold base-rates.
- **Missing per the brief:** distance-based / confidence-weighted negatives and a PU-learning
  comparison are not implemented.

---

## CRITICAL TASK 4 — GNN validation

- **Attribution impossible.** Every "GNN" couples a graph branch with a parallel feed-forward
  `bypass` MLP via a learned gate `alpha` (`gnn.py:137-168`). There is **no graph-only
  condition**, so improvements cannot be attributed to the graph. The brief's required
  tabular-only / graph-only / hybrid decomposition is not implemented.
- **Silent fallback.** `_run_fold()` catches *any* torch exception and silently returns
  HistGradientBoosting predictions still labeled "GAT"/"SAGE" (`gnn.py:370-378`). Reported GNN
  numbers may not be GNN numbers.
- **Transductive check:** message passing uses **train-positive edges only**
  (`gnn.py:174-186, 240`), and in LOSO the held-out host node has no training edges (isolated,
  self-loop only), so neighborhood aggregation does not directly leak the test label — this part
  is acceptable. The real GNN problems are attribution, silent fallback, and global node features
  (Finding C).
- **Alpha interpretation** (`gnn.py:331-336`) uses arbitrary 0.4/0.6 bins; fine as description,
  weak as biology.

---

## CRITICAL TASK 5 — Biological feature validation

- Feature groups exist and are ablatable (`ablation_study._feature_groups`), but the brief's
  **permutation importance, SHAP, and grouped SHAP are not implemented at all.**
- Pair features (`tetra_corr/cub_dist/gc_match/len_ratio`) are computed from real sequences and
  are legitimately available at prediction time (`data_enrichment.py:303-333`) — good.
- PCA/SVD reductions are global (Finding C), so any feature-importance story is biased until
  refit per fold.

---

## CRITICAL TASK 6 — Generalization

| Benchmark | Implemented? | Notes |
|---|---|---|
| LOSO (host species) | Yes | but measures **host** generalization; a phage can be in train+test (disclosed `summary.py:189-193`) |
| LOGO (genus) | Yes | correct counter to genus leakage |
| Unseen-phage MC | Yes (`strain_sim.py`) | reasonable; report CIs |
| Leave-Clade-Out (phylogenetic) | **No** | required by brief; not present |
| Combined unseen host + unseen phage | **No** | the hardest, most informative test; not present |
| Homology-removed re-evaluation | **No** | required to bound Finding E |

---

## CRITICAL TASK 7 — Statistical rigor

- **DeLong** (`statistical_tests.py:90-130`): structurally correct (Mann–Whitney components,
  covariance via `np.cov`), guarded to require identical label vectors
  (`statistical_tests.py:266-268`). Acceptable.
- **McNemar** (`:136-150`): discordant cells correctly enumerated; continuity-corrected. OK.
- **Wilcoxon + BH-FDR** (`:28-84, 230-249`): correct, but BH is applied **only** to Wilcoxon
  p-values; DeLong/McNemar p-values are left uncorrected → inconsistent multiple-comparison
  handling.
- **Bootstrap** (`:156-199`): standard + genus-cluster bootstrap — good, publication-appropriate.
- **Missing (brief):** permutation testing, **calibration analysis (ECE/reliability)**, and
  uncertainty estimation are absent. `metrics()` substitutes AUC=0.5 / PR=0.0 on any failure
  (`evaluation.py:100-159`), silently biasing fold means; PR-AUC via `np.trapz` over the PR
  curve is less accurate than `average_precision_score`.

---

## CRITICAL TASK 8 — Reproducibility

- No `config.yaml`, no data, no `requirements.txt`, no pinned versions, no `DATA.md`, not a git
  repo, Python version unconstrained (see §0). SHA-256 freezing exists (`data_collection.py:476-489`)
  but nothing to hash here.
- Nondeterminism: `RandomForest(n_jobs=-1)` and `XGBoost(n_jobs=-1)` are not bitwise-deterministic;
  no `PYTHONHASHSEED`; torch determinism flags (`use_deterministic_algorithms`, cudnn) not set
  (`gnn.py:198`). The `--thorough` flag silently changes results-affecting hyperparameters
  (`run_pipeline.py:127-132`).
- The Methods/repro card are auto-generated and will read as authoritative even though the
  pipeline can't currently run — a credibility risk.

---

## CRITICAL TASK 9 — Visualization

`plotting.py` is clean (300 DPI, colorblind palette, captions). But it produces a *subset* of
the brief's required figures and several are **methodologically tied to the flawed numbers
above** (e.g. Fig 16 experimental validation = Finding A; Fig 14 leakage audit will look
reassuring while exact-pair-only). Missing: ROC/PR curves, calibration/reliability, infection
network + graph embeddings, SHAP/grouped SHAP, host-range / broad-host-range analyses,
per-genus generalization, LOSO-vs-LOGO-vs-unseen comparison.

---

## CRITICAL TASK 10 — Model improvement

CatBoost / LightGBM / calibrated / stacked ensembles are **not** present (only RF, HistGB,
XGB, GAT, SAGE). This is fine — but no new model should be added until Findings A–E are fixed,
because every metric on the current dataset is contaminated by them. Optimize validity first.

---

## CRITICAL TASK 11 — Code review (software)

- **Pervasive `except Exception` → silent degradation.** Scientific pipelines should fail loudly;
  here failures become AUC=0.5, empty frames, or zero-vectors (`evaluation.py:100-159`,
  `gnn.py:370-378`, many sites). This is the most systemic software risk.
- **`run_experimental_validation`** maps exp-val phages/hosts to matrices and zero-fills unknowns
  (`evaluation.py:330-336`); S. aureus phages are likely absent ⇒ near-zero features ⇒
  meaningless predictions; also `scaler.transform` is called 3× redundantly (`:347-349`).
- **`--frozen-dataset` path bypasses `build_negatives`** (`run_pipeline.py:202-214`); good for
  determinism, but means the frozen CSV must already encode negatives + features consistently, or
  results silently differ from the rebuild path.
- Threaded feature compute preserves order via index map (`data_enrichment.py:108-124`) — correct.

---

## FINAL SCIENTIFIC ASSESSMENT

1. **Are the reported metrics trustworthy?** No — not until Findings A (circular validation), B
   (imputation leak), C (global PCA/SVD), D (name leakage), and E (no homology control) are fixed
   and the numbers re-measured. Expect headline AUCs to **drop** after correction; that is the
   correct outcome.
2. **Are the GNN results biologically meaningful?** Cannot be claimed — no graph-only baseline,
   and silent fallback can relabel gradient boosting as "GNN."
3. **Are the negative-sampling assumptions defensible?** Partially. The open-world caveat is
   stated, but cross-genus easy negatives inflate metrics and no PU/confidence-weighted
   alternative is evaluated. Make `hard_negatives_only` the headline.
4. **Is there evidence of leakage?** Yes — name leakage (D) is structurally present; circular
   validation (A) is certain; imputation leakage (B) is latent-critical; homology leakage (E) is
   uncontrolled.
5. **Is the repository ready for peer review?** No.
6. **What would a Nature Biotechnology reviewer criticize first?** The "independent experimental
   validation" being a relabeled subset of the training source (Finding A), immediately followed
   by the absence of any sequence-homology-aware splitting (Finding E).
7. **Remaining weaknesses to disclose in the manuscript:** open-world negative assumption; LOSO
   measures host- not phage-generalization; no clade-level or combined unseen-host+phage
   benchmark; missing calibration; unverified reproducibility environment.

### Prioritized remediation (do in this order)
1. Decouple experimental validation from the training source, or relabel it honestly (Finding A).
2. Remove/recompute the q75 imputation for pair features; never impute label-correlated
   pair-level columns (Finding B).
3. Move PCA/SVD/name-embedding fitting **inside each training fold** (Finding C).
4. Add sequence-homology dedup (MASH/ANI) + Leave-Clade-Out + combined unseen benchmarks
   (Findings E, Task 6).
5. Make GNN attribution real: add graph-only and tabular-only conditions; make the classical
   fallback **explicit and labeled**, never silent (Finding G).
6. Drop original-name embeddings from the headline model; tighten the leakage threshold
   (Finding D).
7. Add calibration, permutation tests, SHAP/grouped SHAP; correct FDR across *all* tests
   (Tasks 5, 7).
8. Ship `config.yaml`, `requirements.txt` (pinned), `DATA.md`, deterministic flags, and a git
   repo (Task 8).

*Audit performed by static reading of all 16 modules; empirical magnitudes pending data + a
Python ≥3.10 scientific environment.*
