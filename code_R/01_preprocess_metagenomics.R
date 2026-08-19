# 01_preprocess_metagenomics.R
# MARS_IBS_2020 (Borenstein curated) -> per-subject CLR matrix
# Input : data_raw/02_Metagenomics/MARS_IBS_2020/{species.tsv,metadata.tsv}
# Output: data_processed/02_Metagenomics/clr_species_n75.rds

library(data.table)

sp <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/species.tsv",
            sep = "\t", data.table = FALSE)
mt <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/metadata.tsv",
            sep = "\t", data.table = FALSE)
stopifnot(nrow(sp) == 444, ncol(sp) == 11514, nrow(mt) == 444)

# 1. collapse to per-subject arithmetic mean of relative abundances
idx  <- match(sp$Sample, mt$Sample); stopifnot(!anyNA(idx))
subj <- as.character(mt$Subject[idx])
M    <- as.matrix(sp[, -1, drop = FALSE])
agg  <- rowsum(M, subj) / as.vector(table(subj)[rownames(rowsum(M, subj))])
stopifnot(nrow(agg) == 75)

# 2. prevalence filter: present in >= 10% of participants (>= 8 of 75)
keep <- colSums(agg > 0) >= ceiling(0.10 * nrow(agg))
agg  <- agg[, keep, drop = FALSE]
stopifnot(ncol(agg) == 4795)

# 3. CLR, zeros replaced by half the smallest positive value
pc <- min(agg[agg > 0]) / 2
L  <- log(ifelse(agg == 0, pc, agg))
X  <- L - rowMeans(L)
stopifnot(max(abs(rowMeans(X))) < 1e-12, all(is.finite(X)))

dir.create("data_processed/02_Metagenomics", recursive = TRUE,
           showWarnings = FALSE)
saveRDS(X, "data_processed/02_Metagenomics/clr_species_n75.rds")
cat("saved:", dim(X), "\n")
