# Lineage-level *Streptococcus* enrichment in IBS-D

Analysis code for a secondary analysis of the MARS_IBS_2020 cohort
(Mars et al., Cell 2020) as curated by the Borenstein laboratory.

## Reproducing the analysis

Scripts are numbered and run in order from the project root:

| Script | Purpose | Self-check |
|---|---|---|
| `01_preprocess_metagenomics.R` | raw table to per-subject CLR matrix | reproduces the stored matrix to 3.6e-15 |
| `02_final_model.R` | final model, Group + Age, n = 74 | reproduces MaAsLin2 coefficients to 6.0e-15 |
| `03_verify_BH_scope.R` | establishes the correction family empirically | 14,385 to 3.0e-15; 9,590 rejected |
| `04_reviewer_response.R` | confirmatory contrasts (A1, A3, A1b) | — |
| `05_external_metadata.R` | detection, depth, PPI, batch | — |

Each script asserts its expected dimensions with `stopifnot()` and fails
loudly rather than proceeding on mismatched input.

## Data not included here

| File | Source |
|---|---|
| `species.tsv`, `metadata.tsv` | github.com/borenstein-lab/microbiome-metabolome-curated-data (`MARS_IBS_2020`) |
| `mmc1.xlsx` (Table S1) | doi.org/10.1016/j.cell.2020.08.007, Supplemental Information |
| `ENA_PRJEB37924_runs.tsv` | ENA Portal API; URL at the foot of `05_external_metadata.R` |

Place them under `data_raw/02_Metagenomics/` as referenced in the scripts.
The processed CLR matrix is included so that `02`-`04` run without the
24 MB raw table.

## Notes on the analysis

- Row alignment is name-based throughout (`match()` on participant ID).
  Positional and logical indexing were the source of an earlier error in
  which the wrong participant was silently excluded.
- The primary model corrects across all 14,385 coefficients MaAsLin2
  returns; the confirmatory analyses correct within each contrast
  (4,795 features). Both scopes are stated where they apply.
- Sequencing batch is unbalanced across groups but cannot be adjusted for:
  IBS-D has almost no within-group variation in it. See the diagnostic in
  `05_external_metadata.R`.

## Environment

`code_R/sessionInfo.txt`. R 4.5.3, Bioconductor 3.21.

## Citation

Wadi WF, Razila S. *Manuscript under review.*

Wisam Faisal Wadi — ORCID [0009-0009-8814-4060](https://orcid.org/0009-0009-8814-4060)
Advanced Medical and Dental Institute, Universiti Sains Malaysia.
