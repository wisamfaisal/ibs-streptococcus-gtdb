# Corrections and verification log

Numerical corrections made to the analysis and its documentation during
pre-submission review, 19 August 2026. Every figure below was read from a
project file or from an external source retrieved in the same session;
outputs are stored under `results/05_Reviewer_Response/`.

---

## 1. Multiple-testing scope: 9,590 was wrong, 14,385 is correct

The working documentation recorded Benjamini-Hochberg correction across
4,795 features x 2 group contrasts = 9,590 tests. MaAsLin2 in fact corrects
across every coefficient it returns, including the age term: 4,795 x 3 =
14,385. Verified in `code_R/03_verify_BH_scope.R`: the 14,385 scope
reproduces the stored q-values to 3.0e-15, whereas the 9,590 scope differs
by 0.114 and is rejected. Fifty coefficients pass q < 0.05 -- 43 in IBS-D,
6 in IBS-C, and one age term (*Collinsella* sp900549455, beta = -0.69,
q = 0.044). The 49 reported in the manuscript are the group-level results.

## 2. Sample-count discrepancy is only partly explained by the stated filter

The original study reports 474 stool samples from 77 participants; the
curated tables contain 444 from 75. Retrieval of ENA PRJEB37924 gives 604
runs: 532 single-end stool libraries and 72 paired-end biopsy libraries,
four of the stool entries being sequencing controls. Of the remaining 528,
430 of the 444 curated samples matched by participant and timepoint. Ninety-
eight archived stool samples are absent from the curated tables, of which
only 16 fall below the repository's stated 50,000-read threshold; the
mechanism excluding the other 82 is not documented. Four curated samples
have no archive counterpart (52.T.6, 50.T.6, 63.T.6, 59.T.2).

## 3. Instrument field in ENA conflicts with the source publication

ENA records `Illumina HiSeq 4000` for all 604 runs. The source publication
reports HiSeq 2500 Rapid Mode (100 bp single-end) or NextSeq (150 bp single-
end). Computed read length (base_count / read_count) is 101 nt for 484 stool
samples, 126 nt for 48, and 301 nt for the biopsies. A single instrument
value across libraries of differing length and layout is a deposition
default rather than a per-run description; the published methods stand.

## 4. Sequencing batch is unbalanced but cannot be adjusted for

Batch assignment, recorded in the original study's supplementary metadata,
differs across groups (Kruskal-Wallis p = 0.006). It cannot be included as
a covariate: the within-group standard deviation of the batch-1 fraction in
IBS-D is 0.115, and 0.082 once PPI users are excluded. Adding it removes all
129 group-level results, and the batch term alone accounts for 774 features
-- a pattern of non-identifiability rather than of a confounder being
corrected. Restricting instead to batch-1 participants (n = 39) leaves all
18 lineages with confidence intervals excluding zero while the matrix-wide
count falls to one, so the signature is not a batch artefact.

## 5. Proton pump inhibitor use does not explain the association

Reported by 3 of 23 IBS-C participants, 3 of 29 IBS-D, and none of the 24
controls. Equal between subtypes, so it cannot account for the difference
between them. Excluding the six users (n = 68) leaves all 18 lineages
significant. Three of the six also reported antibiotic use, so the two
sensitivity analyses overlap.

## 6. Taxonomic resolution of the comparison reference

The documentation recorded 26 *Streptococcus* species in the reference used
by Vich Vila et al. Counting directly from their Table S2 gives 22 species,
among 302 species-level features in total, five of the 22 being unnamed or
environmental designations.

---

## Methodological note

Two documentation lines were briefly overwritten during this session because
line numbers were used as edit targets without re-reading the file. The file
was restored from a copy and edited again by text match, asserting exactly
one occurrence before writing. This is the same failure mode recorded
earlier for participant exclusion, where a logical mask silently removed the
wrong participant: address rows and strings by name, never by position.
