# Original-source comparator recovery: public scientific summary

Audit date: 15 October 2026. This is a retrospective public-resource feasibility
audit, not new outcome-blind validation or a completed predictor rerun.

The exact all-100 continuous-score comparison remains unavailable. The released
decision table contains 300 selected constituent records: 233 numeric scores
and 67 absent values, all at the generalist `4_rescue` stage. Exactly 64 hosts
have three numeric scores. The rescue procedure does not establish that a
host-specific numeric probability existed for every selection.

The current upstream main revision matches the pinned source commit
`99656d8c020fcd1dbbc30a02ee10274a07655562`. The audit inspected the released tree,
pointwise CV script, external feature coverage and every one of the 96 CSV
members in the released prediction RAR (1,918,900 rows). Those predictions cover
402 development host names, not authentic external-host predictions. The sole
shared name is AN03, an identifier collision involving different genomes.

No fitted pointwise classifier/preprocessing object or runnable end-to-end
four-stage recommender was released. The pointwise script's preflight stopped
before fitting on an unavailable original training-matrix file. External capsule
fields and 99 query UMAP coordinates/the fitted projection were also absent.
Replacing those resources or assigning generalist coverage as a missing
probability would create a reconstruction, not recover original source scores.

The previous matched 64-host comparison and missing-score identification bounds
remain the available evidence. No missing value was imputed, model fitted,
downloaded pickle executed or upstream source changed. No email or new issue was
sent during this audit. Personal commit metadata, contact drafts and original
private machine-path details are excluded from the public summary.

The saved records include score coverage by stage, the biological record IDs
with missing scores, prediction-archive coverage and the retrospective decision
protocol. Full reproduction of the preflight requires a pristine pinned source
checkout with the RAR; the compact CSV subset alone is not that checkout.

Sources: [authors' pinned repository](https://github.com/mdmparis/coli_phage_interactions_2023/tree/99656d8c020fcd1dbbc30a02ee10274a07655562)
and [Gaborieau et al.](https://doi.org/10.1038/s41564-024-01832-5).
