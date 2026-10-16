# BASEL–ECOR uncertainty specification before outcome evaluation

Recorded after genome-only overlap auditing, before the new predictions were
joined to experimental labels. At the historical 0.05 canonical Mash-form
threshold, all 25 ECOR genomes belong to one connected host cluster; the 52
phages form 31 clusters. Therefore a two-way bootstrap at that threshold cannot
estimate between-host-cluster generalization uncertainty. Report that primary
interval as unavailable, rather than treating repeated copies of one cluster
as independent information.

Retain the frozen point estimate and every eligible pair regardless of result.
As explicitly secondary, conditional uncertainty diagnostics, report the same
2,000-replicate two-way bootstrap with host clusters at both 0.001 and 0.01,
keeping phage clusters at 0.05, plus an entity-identity two-way bootstrap. Report
all three definitions, group counts, skipped/degenerate replicates, and any
undefined interval. These narrower host distances distinguish strains within
the single genus; they do not validate generalization to new species. No
distance definition is to be selected based on whether its interval excludes
chance. Sketches are canonical, contig-aware, k=21 with 2,048 hashes and seed 42.

All 25 external hosts are within 0.05 of a training host reference. Consequently
this test cannot establish dual-axis sequence-cold transfer at that cutoff.
Report source-held-out status, exact contig-set overlap, and sequence-nearness
separately. Absence of exact contig-set equality does not prove strain novelty,
particularly when assemblies differ in fragmentation.
