# Additional-bank and canonical-sketch analysis plan

Recorded 2026-10-05 before running the analyses below. This is a retrospective
extension after previous PrecisionPhage development, not a preregistration or
a new blinded test. Preserve original frozen predictions and input files.

## Canonical genome grouping

Use sourmash 4.9.3 DNA MinHash (k=21, seed=42), processing each FASTA record
separately, skipping ambiguous windows, with 256 and 2,048 hashes. Confirm
reverse-complement invariance. Transform estimated Jaccard similarity with
the Mash distance formula and form connected components at distance <=0.05.
These remain sketch distances, not measured ANI. Compare cluster membership
with the historical forward-strand implementation and record FASTA checksums.

For the primary sensitivity, retain exactly the historical 398 test pairs and
five blocks. Use the 2,048-hash canonical clusters to define shared-entity,
host-excluded, phage-excluded, and both-excluded training pools. Exclude test
rows in every regime. Match positive and negative training counts within each
block across the four regimes. Use the inherited 400-tree configuration and
seeds 20261005, 20261006, 20261007. Save fold results, predictions, overlap
checks and non-estimable folds. Do not select a favourable sketch size or seed.

## Replication of the selection workflow across additional banks

Use the released GenoPHI Klebsiella1, Klebsiella2, Pseudomonas and Vibrio
interaction matrices with their published host phylogenies where exact
identifiers resolve. Exclude the Escherichia development matrix from claims
of additional-bank replication. Document source URLs, licenses, checksums,
duplicate rows, missing labels and tree-tip coverage before modeling.

Freeze the nearest-neighbour model family and k grid (5,10,20,40), selecting
k within three inner grouped folds of each of five outer grouped folds. Fit
separately within each bank. Primary selection is the three highest raw
distance-weighted neighbour scores. Metadata constraints from the Picard bank
cannot be assumed for other banks; report this as an unconstrained adaptation,
not validation of the complete original workflow or a frozen-model transfer.
Stable identifier ordering breaks ties. Evaluate all hosts, with missing
assays represented by lower/upper coverage bounds.

Compare with training-only phage prevalence and a stronger greedy observed
coverage triplet. Use host identity groups and tree patristic connected
components at 1e-4, 1e-3, 1e-2; branch lengths are source-specific and are not
directly interchangeable with the Picard distance scale. Skip settings with
fewer than five groups or fewer than two inner training groups. Report every
eligible bank and threshold, five-fold coverage and paired group-bootstrap
95% intervals (2,000 resamples, seed 20261005), plus AUROC/AP over measured
pairs as descriptive discrimination metrics. Group bootstrap conditions on
the bank and fitted cross-validation predictions and does not include all
training uncertainty. Do not pool banks as independent replicates of a single
population or claim mixture efficacy from individual constituent activity.

If published phylogenies are unavailable or identifiers cannot be resolved
exactly, report the concrete missing input and do not invent genome mappings.
All work uses local CPU resources; no paid compute is needed for this plan.
