# Source provenance and third-party terms

PrecisionPhage code is supplied under the accompanying MIT LICENSE. This does
not relicense third-party datasets or software.

- Gaborieau et al. (2024), Nature Microbiology 9:2847-2861,
  DOI https://doi.org/10.1038/s41564-024-01832-5. Compact original CSV inputs are
  unchanged from https://github.com/mdmparis/coli_phage_interactions_2023 at
  commit 99656d8c020fcd1dbbc30a02ee10274a07655562. Its archived code/data record
  https://zenodo.org/records/13831957 identifies Creative Commons Attribution
  4.0 (https://creativecommons.org/licenses/by/4.0/); checked through the primary
  record API during packaging. Derived PrecisionPhage analyses are identified
  separately by their paths. No article PDF is redistributed.
- VHIP: Bastien et al. (2024), DOI https://doi.org/10.1371/journal.pcbi.1011649.
  Upstream analysis commit 5c87aaf7164f78c11f598e0bb56681932ce8bbad and tool
  commit b62a7b6ba08a056b6e04b25ef4dd3a355ebd179d. The tool's MIT license is
  included under reproducibility/third_party_licenses. The input matrix is
  retained with exact hashes; primary source terms continue to apply.
- Goller et al. (2021), DOI https://doi.org/10.1038/s41467-021-27037-6,
  is the experimental source underlying StaphStudy. See the VHIP provenance
  and the source-held-out protocol for its use here.
- Noonan et al. (2026), DOI https://doi.org/10.1038/s41564-026-02482-5,
  provided the BASEL-ECOR labels. GenoPHI commit
  f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f pins that CSV. Its license text is
  provided separately where available; no GenoPHI fitted-model collection is
  distributed. Additional-bank archives are obtained with the optional fetch
  command and retain their source terms.
- BASEL phage accession table: https://doi.org/10.1371/journal.pbio.3003063.s002.
  ECOR assembly accession table: https://doi.org/10.1128/MRA.01133-18. The
  XLSX accession table and derived accession mapping accompany the frozen
  analysis. The original ECOR article XML is not redistributed: the existing
  preparation script retrieves its primary Europe PMC URL when required, and
  its historical checksum remains in genome_accession_mapping.json. This
  avoids redistributing an entire article and unneeded correspondence details.
- ENA accession OZ035781.1 establishes the CLB_P2/55989_P2 identity through
  exact sequence correspondence, including a documented 60-base terminal
  prefix duplication in the local source assembly. The FASTAs themselves are
  fetched only when the source-score audit is rerun.
- PHIST: https://github.com/refresh-bio/PHIST, pinned local source commit
  849efbe1894e321c1ff06b13f7f1aca377bfddb5. Its GPL license is included as an
  attribution record. PHIST source/build artifacts are not bundled here; the
  saved aligned predictions permit the reported metric audit without a build.

The original root README/DATA prose is not used as the current scientific
description. Historical protocols/logs can contain superseded interpretations.
The revision README and revised manuscript define the current claims.
