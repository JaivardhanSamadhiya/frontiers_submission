# Additional-cohort and authentic-source audit, 17 October 2026

## What was checked

The public source repository HEAD was rechecked with `git ls-remote` on 17 October 2026. It remains `99656d8c020fcd1dbbc30a02ee10274a07655562`, identical to the locally audited Gaborieau source. No newer public source commit supplies the 67 rescue-selection probabilities. The earlier archived-data audit is retained. This does not establish that the authors never generated those scores. The author reports sending the missing-score request and attachment on 16 October; no reply or new private scores have been supplied to this analysis. No duplicate email was sent.

Source: https://github.com/mdmparis/coli_phage_interactions_2023

## Eligibility for another test of the repaired fixed-bank scorer

The existing Picard-trained scorer contains fitted models for specific phages. Another Escherichia dataset is not automatically compatible: the tested phages need exact or independently verified identity matches to those fitted bank members. A qualifying cohort also needs all required host covariates and genomic-distance inputs, measured host–cocktail identities and outcomes, and a provenance record establishing that those outcomes did not guide this scorer's repair or selection. Existing ColoColi observations and BASEL-ECOR outcomes have already been examined and cannot become untouched by changing the split or random seed.

A new bank can support a separate, transparently retrained fixed-bank experiment, but cannot serve as frozen external validation of the existing 96-phage bank model without demonstrated compatibility. Broad genome downloading does not resolve this distinction.

## Primary resources examined

1. **Gaborieau et al. (2024), Picard/ColoColi.** This supplies the same 100-host cocktail cohort already analyzed. The currently visible source commit is unchanged. The newly analyzed doublets are useful retrospective decisions within this resource, not a second independent cohort.
   - https://doi.org/10.1038/s41564-024-01832-5
   - https://github.com/mdmparis/coli_phage_interactions_2023/tree/99656d8c020fcd1dbbc30a02ee10274a07655562
2. **Noonan et al. (2026), GenoPHI experimental validation.** Public experimental matrices and genome accessions exist. The local study has already used the BASEL-ECOR matrix, and BASEL phages are not thereby established as the Picard model's fixed phages. This resource cannot be described as another untouched cocktail test in this project. It remains relevant to the sequence-model external benchmark and comparator work already documented.
   - https://doi.org/10.1038/s41564-026-02482-5
   - https://doi.org/10.5281/zenodo.21696332
   - https://github.com/Noonanav/GenoPHI
3. **Predictive phage therapy for Escherichia coli urinary tract infections (2024).** The study reports 31 phages, 314 hosts, a later ten-host evaluation, and experimentally tested cocktails in artificial/canine urine. This is a potential separate-bank replication resource. However, no verified mapping to the current bank or ready compatible inputs for its frozen models was established. The phenotype definition and environment also differ. Training a fresh model on this resource would answer a different question and requires a separately locked protocol and input audit before inspecting its held-out numerical results.
   - https://doi.org/10.1073/pnas.2313574121
   - https://www.pnas.org/doi/pdf/10.1073/pnas.2313574121
   - Genome resource reported by the authors: https://enterobase.warwick.ac.uk/species/ecoli/search_strains?query=workspace:97040
4. **Penke et al. (2026), automated therapeutic-cocktail screening.** This reports a large UPEC screening program and cocktail assays, so it is relevant scientifically. The paper describes a proprietary cocktail predictor and a different phage bank. The inspected public material did not establish exact phage mapping, the complete compatible feature/assay package, or a second test of this fitted scorer. It was not treated as validation merely because a large screening count was reported.
   - https://doi.org/10.1038/s41467-026-68684-x
   - https://escholarship.org/uc/item/2kd1z89x

## Result and concrete next requirement

This scoped search identified plausible separate-bank research resources but **did not establish an eligible, untouched outcome cohort for the repaired Picard scorer**. This is a compatibility/access limit, not a result of model performance. No new-cohort performance claim is made and no numerical outcomes from a putative untouched cohort were downloaded to select an attractive test.

Before such a test, obtain a host–cocktail manifest with verified fixed-bank phage identities, the required host input data, assay definitions, and an outcomes file kept separate from prediction generation. Freeze scorer/artifact hashes, inclusion rules, ranking formula, tie handling, baseline definitions and host/composition clustering. Record predictions and hashes before joining labels. If those resources cannot be obtained, retain the retrospective limitation rather than reclassifying existing observations as new validation.

## Comparator execution

`comparator_resume_genophi.py` resumes the earlier fully budgeted, **adapted** GenoPHI comparison in `data/results_v2/frontiers_strengthening_2026_10_16/comparator/official_full_budget_resumed`. It preserves original files, budgets and the disclosed negative-MCC retention repair; the helper changes output location and worker dispatch only. Five outer folds and the existing independent audit must finish before reporting comparative performance. This remains a Pseudomonas known-bank, held-out-host comparison; it is not the missing Gaborieau comparison or a Picard-to-new-bank experiment.
