# Scientific snapshot v2026.10.15.1

## Archive-portability patch

The patch records the exact staged Git-blob SHA-256 of each inventoried file,
in addition to its working-tree hash and the existing explicit text-EOL policy.
Git-generated archives can therefore be checked without treating unlisted
formats as normalized text. The verifier accepts only recorded hashes; it
does not allow arbitrary file transformations. CI now also runs the public
payload and saved-prediction verifiers. LF normalization is explicitly recorded
for the known text configuration names `.gitignore`, `.gitattributes` and
`.gitmodules`, covering Windows-generated archives without treating binary
files as text. The original `v2026.10.15` tag remains
unchanged. Scientific inputs, predictions, bootstrap results and claims are
unchanged by this packaging correction.

This dated release preserves the v2.0.1 public repository history and original
frozen predictions while adding retrospective scientific audits. It is not a
new prospective validation, clinical product or independently blinded trial.

## Additions and corrections

- Explicit AP prevalence references, matched evaluation cohorts, corrected
  individual fold metrics and coverage-bias tables.
- Training-only degree controls; phage-only, host-only, pair-only and full
  feature ablations; canonical sketch/cutoff sensitivities; a cache-limited
  unsuccessful Nahant sequence-resolution audit.
- Frozen negative BASEL–ECOR transfer results and retrospective target-label
  recalibration learning curves, kept separate from raw frozen transfer.
- Trivial cocktail controls, one-source 64-complete-host comparison, missing-score
  bounds, assay-triage summaries and a batch-invariant retrospective scoring
  repair. No missing source rescue probability is imputed.
- Crossed host-cluster × selected-composition bootstrap sensitivities with
  saved weights and paired draws, plus reproduction of all 1,000 whole-triplet
  shuffles. The shuffle preserves menu repetition but remains a diagnostic,
  not an exchangeability-justified significance test.
- A paired two-measured-options dependence sensitivity. The policy-gain interval
  includes zero; choosing new cocktail constituents is not validated.
- A documented attempt to recover the original-source all-100 predictor. Its
  released resources do not support a faithful complete rerun. Public records
  retain the failed feasibility gate, not private contact drafts.
- Complete frozen training-genome accession/order/hash provenance and an
  opt-in, non-overwriting recovery helper. Large sequences remain excluded.
- Exact public-payload checksums, actual revision-audit package versions,
  study-only reproducibility code and citation/archive metadata.

## Preserved limitations

Original tree-model selection rationale and exact historical tuning settings
cannot be fully reconstructed. Earlier tuning logs use a different feature
schema. Broad source transfer, portable calibrated probabilities, superiority
to unrun released tools, meaningful within-host selection improvement,
resistance suppression and clinical efficacy are not established.

The historical cocktail AUROC 0.807 and AP 0.951 are unchanged. The host-only
AUROC interval 0.693–0.906 is conditional on the tested menu; the later crossed
sensitivity is wider, 0.543–0.965. Neither supports transport to a new bank or
prospective selection procedure. Source decision scores cover 64 complete hosts,
not all 100; the exact full comparison remains unavailable.

Writing/submission automation is removed from the new tip without rewriting
Git history. Historical scientific code and outputs remain labelled as history
where required for auditability. No manuscript or private peer-review material
belongs to this scientific payload.
