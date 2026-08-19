# 02_final_model.R
# Final model: MaAsLin2-equivalent OLS, Group + Age, n = 74
# Input : data_processed/02_Metagenomics/clr_species_n75.rds
#         data_raw/02_Metagenomics/MARS_IBS_2020/metadata.tsv
# Output: results/05_Reviewer_Response/model_final_n74.rds
# Note  : participant 10007545 (sole male IBS-C) excluded post hoc,
#         see MODEL_SELECTION_README. BH applied jointly across all
#         features and both group contrasts (9,590 tests).

library(data.table)

X75 <- readRDS("data_processed/02_Metagenomics/clr_species_n75.rds")
mt  <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/metadata.tsv",
             sep = "\t", data.table = FALSE)

md <- unique(mt[, c("Subject", "Study.Group", "Age")])
md <- md[match(rownames(X75), as.character(md$Subject)), ]
rownames(md) <- as.character(md$Subject)
stopifnot(identical(rownames(X75), rownames(md)), !anyNA(md))

keep <- rownames(X75) != "10007545"
X <- X75[keep, , drop = FALSE]; M <- md[keep, ]
stopifnot(nrow(X) == 74, identical(rownames(X), rownames(M)))

M$Group <- relevel(factor(M$Study.Group), ref = "H")
M$Age_z <- as.numeric(scale(M$Age))
M$Y <- X

fit <- lm(Y ~ Group + Age_z, data = M)
mm  <- model.matrix(fit); B <- coef(fit)
iv  <- solve(crossprod(mm)); dfr <- nrow(M) - ncol(mm)
s2  <- colSums(residuals(fit)^2) / dfr

pull <- function(k, lab) {
  se <- sqrt(iv[k, k] * s2); b <- B[k, ]; tt <- b / se
  data.frame(feature = colnames(X), value = lab, coef = b, stderr = se,
             ci_lo = b - qt(.975, dfr) * se, ci_hi = b + qt(.975, dfr) * se,
             pval = 2 * pt(-abs(tt), dfr), row.names = NULL)
}
res <- rbind(pull("GroupD", "D"), pull("GroupC", "C"))
stopifnot(nrow(res) == 9590)
res$qval <- p.adjust(res$pval, "BH")

dir.create("results/05_Reviewer_Response", recursive = TRUE,
           showWarnings = FALSE)
saveRDS(res, "results/05_Reviewer_Response/model_final_n74.rds")
cat("significant (q<0.05):", sum(res$qval < 0.05), "\n")
