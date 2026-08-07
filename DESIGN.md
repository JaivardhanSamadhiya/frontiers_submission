# PrecisionPhage v2 — Architecture & Design Decisions

A ground-up, reproducibility-first, leakage-free phage–host interaction framework built around
**real experimental cross-infection data**, a leakage-safe bipartite GNN + Edge-MLP decoder,
rigorous generalization testing, cocktail optimization, and temporal resistance modeling.

This document records, for each component, the candidate approaches considered and the one
chosen, with rationale grounded in the data we actually have and the current literature
(CHERRY, HostG, PHPGAT, iPHoP; strain-level review PMC12936788; ViroBench temporal-split
findings). Each implemented module is reviewed against this document.

---

## 0. The decisive reframe — use the real labels

`data/raw/VirusHostInter.csv` contains **8,849 experimentally assayed pairs with true labels**:

| Study | Rows | Nature |
|---|---|---|
| NahantCollection | 5,208 | Vibrio dense cross-infection matrix (Kauffman et al., *Nature* 2018) |
| NCBI_HR | 2,588 | NCBI curated host-range |
| StaphStudy | 1,053 | *Staphylococcus* infection matrix |

Totals: **2,770 Inf / 6,079 NoInf**, 2,331 phages, 392 hosts; precomputed pairwise features
(k3dist, k6dist, GCdiff, Homology) present for **both** classes.

**Decision:** the supervised model is trained and evaluated **only on real experimentally-labeled
pairs**. Constructed negatives are abandoned as the primary signal (they caused the v1 collapse:
hard-negative AUC 0.24, i.e. anti-predictive, because within-genus "negatives" were
false-negative-contaminated). This single decision removes the project's biggest leakage/validity
threat. Constructed negatives may appear later *only* as an explicitly-labeled robustness ablation
under PU-learning assumptions, never in headline metrics.

---

## 1. Data layer

- **Considered:** (a) re-merge all 5 heterogeneous sources as in v1; (b) center on the dense
  experimental matrices; (c) supplement with iPHoP/VHDB.
- **Chosen: (b)** with optional (c) as external transfer test. Dense matrices give real negatives,
  enable strain-level prediction, and are exactly what cocktail optimization needs (a phage×host
  coverage matrix). The three studies are kept as a `study` column so we can do
  **cross-study (cross-laboratory) transfer** — the strongest generalization test per the
  literature.
- **Genome linking:** each phage/host name is resolved to a FASTA via an explicit, cached
  name→accession map (no fuzzy slug guessing in the hot path). Pairs lacking a resolvable genome
  for *either* partner are excluded from sequence-based models and reported, not zero-imputed.

## 2. Leakage-safe splitting (the core scientific contribution)

Every reviewer concern in the field is about homology leakage. We implement **grouped, similarity-aware
splits** and fit *all* preprocessing inside the training fold.

- **Similarity clustering:** cluster phages by genome similarity (MinHash/MASH-style sketch +
  greedy single-linkage at a configurable ANI/Jaccard threshold) and hosts likewise. `sourmash`/
  Mash-style sketching in pure-Python/Bio if no external binary.
- **Evaluation regimes (all reported):**
  1. **LOSO** — leave-one-host-species-out (host-species generalization).
  2. **LOGO** — leave-one-host-genus-out.
  3. **Leave-one-phage-cluster-out** — unseen phage (the test v1 lacked).
  4. **Combined unseen-phage × unseen-host** — hardest; the headline novelty claim.
  5. **Cross-study transfer** — train on two studies, test on the third (domain shift).
- **Invariant:** the message-passing graph, scalers, PCA, and k-mer vocabularies are built from
  **training rows only**, enforced by a single `Fold` object that exposes train indices and a
  `fit_transform`/`transform` contract. A unit test asserts no test node/edge influences training
  tensors.

## 3. Genomic features

- **Phage/host node features:** k-mer spectra (k=4 canonical, the CHERRY-proven representation),
  codon usage, GC, genome length, dinucleotide bias. Reduced by **fold-internal** PCA.
- **Edge (pairwise) features → Edge-MLP:** features computable for *any* candidate pair so they
  never leak: k-mer distance (d2*-style), GC difference, oligonucleotide-frequency correlation,
  and alignment homology (BLASTN coverage/identity if `blastn` available, else k-mer proxy).
  These feed the **decoder MLP**, which is the user-requested "Edge MLP" integration point.
- **No name embeddings** in the headline model (v1 showed they were net noise and a leakage risk).

## 4. Model — bipartite GNN encoder + Edge-MLP decoder

- **Considered:** GCN (transductive, CHERRY/HostG), GAT, **inductive GraphSAGE**, plain MLP.
- **Chosen: inductive GraphSAGE encoder + Edge-MLP decoder.**
  - Encoder message-passes over **training-positive** phage–host edges to produce node embeddings;
    inductive so that an unseen host/phage (no training edges) still gets an embedding from its
    node features. This directly fixes v1's failure mode (in LOSO the held-out node was isolated,
    `alpha` froze at 0.5, the graph added nothing).
  - **Decoder:** `MLP([z_phage ‖ z_host ‖ edge_features]) → P(infection)`. Because the decoder
    always sees genomic pairwise features, predictions are informative even for graph-isolated
    nodes — the graph *adds* relational signal rather than being a single point of failure.
  - **Baselines for honest attribution (all run):** (i) Edge-MLP on features only (no graph),
    (ii) GraphSAGE embeddings only (no edge features), (iii) full hybrid, plus gradient-boosting
    (XGBoost/HistGB) on the flat feature vector. The graph's value = hybrid − features-only,
    measured under every split regime.
- **Determinism:** fixed seeds, `torch.use_deterministic_algorithms(True)`, capped threads,
  pinned versions.

## 5. Negatives

- **Primary:** real experimental `NoInf`.
- **Robustness ablation only:** PU-learning (elkan-noto style) and phylogeny-distance negatives,
  reported separately with explicit assumptions; never mixed into the headline number.

## 6. Evaluation & statistics

- Metrics with bootstrap CIs (genus-clustered), **calibration (reliability curve + ECE)**,
  DeLong, permutation tests, McNemar. No silent 0.5 substitution — degenerate folds are recorded
  and excluded transparently. Per-regime tables + a single results manifest with SHA-256.

## 7. Cocktail optimization

- Input: predicted P(infection) matrix over a target bacterial population (strains).
- **Algorithm:** weighted maximum-coverage / set-cover — greedy (with 1−1/e guarantee) and exact
  ILP (PuLP/CBC) for small instances; objective = strains covered above a calibrated probability
  threshold, with diversity/breadth regularization.
- Output: minimal robust cocktail + coverage curve, validated against held-out infection matrices.

## 8. Temporal resistance modeling (the *Nature* differentiator)

- **Considered:** static breadth only; pure ODE eco-evolution; GNN-coupled eco-evolution.
- **Chosen: GNN-coupled eco-evolutionary simulation.**
  - Population dynamics ODE/stochastic model: susceptible + resistant bacterial subpopulations and
    phage populations (resistance-cost, mutation-rate, burst-size parameters from literature).
  - The GNN predicts cross-resistance structure: when a strain evolves resistance to phage A, which
    other phages still infect it (orthogonal receptors ⇒ low predicted cross-resistance).
  - Deliverables: **time-to-resistance curves**, and a **resistance-robust cocktail design** that
    optimizes long-horizon coverage, not just day-0 breadth. This converts a static predictor into
    a therapeutic-design tool — the headline application novelty.

## 9. Repository / reproducibility

```
full_pipeline/
  pyproject.toml            # pinned deps, console entrypoints
  configs/default.yaml      # all knobs; no logic in config
  src/precisionphage/
    data/  features/  splits/  models/  eval/  cocktail/  temporal/  viz/  utils/
  experiments/              # thin runnable drivers
  tests/                    # unit + leakage invariants
  DESIGN.md  (this file)
```

- Deterministic, seeded, version-pinned; frozen dataset + SHA-256; every figure regenerable from
  one command; leakage invariants enforced by tests in CI.

---

## Self-review checkpoints (performed after each module)
1. Does any test-fold information touch training tensors? (must be: no)
2. Is every reported feature computable at prediction time for an unseen pair? (must be: yes)
3. Are negatives real or clearly-labeled assumptions? (must be: real for headline)
4. Does the graph's contribution survive the features-only baseline under the hardest split?
5. Is the result reproducible from a clean checkout with one command?
