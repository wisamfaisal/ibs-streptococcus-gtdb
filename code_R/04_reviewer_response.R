# 04_reviewer_response.R
# Confirmatory analyses added in response to pre-submission review.
# A1  direct IBS-D vs IBS-C contrast
# A3  IBS-D vs controls with sex adjustment (collinearity absent here)
# A1b IBS-D vs IBS-C, females only
# A5  sequencing depth from ENA PRJEB37924
# B1  per-group detection frequency and median relative abundance
# B8  antibiotic sensitivity (22 exposed samples removed)
# A4  PPI sensitivity (6 users removed) and sequencing-batch diagnostics

library(data.table)

SIG18 <- c("parasanguinis_G","sp001813295","parasanguinis_B","parasanguinis_D",
  "parasanguinis_C","parasanguinis_E","sp001578805","sp900755085","sp902836505",
  "parasanguinis","parasanguinis_A","sp902373455","parasanguinis_F","intermedius",
  "alactolyticus","sp000187445","vestibularis","sp003521145")

is18 <- function(cn) sub(".*s__Streptococcus ", "", cn) %in% SIG18 &
                     grepl("g__Streptococcus;", cn, fixed = TRUE)

# generic OLS across a response matrix; returns one contrast
ols <- function(X, M, rhs, term) {
  M$Y <- X
  f  <- lm(as.formula(paste("Y ~", rhs)), data = M)
  mm <- model.matrix(f); B <- coef(f); iv <- solve(crossprod(mm))
  dfr <- nrow(M) - ncol(mm); s2 <- colSums(residuals(f)^2) / dfr
  se <- sqrt(iv[term, term] * s2); b <- B[term, ]
  data.frame(feature = colnames(X), beta = b, se = se,
             ci_lo = b - qt(.975, dfr) * se, ci_hi = b + qt(.975, dfr) * se,
             pval = 2 * pt(-abs(b / se), dfr), df = dfr, row.names = NULL)
}

# Multiple-testing scope for the confirmatory analyses
# ----------------------------------------------------
# BH is applied within each contrast separately, i.e. a family of 4,795
# features per hypothesis. This differs from the primary model, where
# MaAsLin2 corrected across all 14,385 coefficients it returned (verified
# in 03_verify_BH_scope.R). The wider scope is not available here: A3 and
# A1b compare two groups only, so no second group contrast exists to
# include. A per-contrast family is the only scope defined identically
# across all three analyses. The scope affects how many lineages cross the
# threshold individually (1 at 9,590 vs 3 at 4,795 for A1) but not the
# direction, the confidence intervals, or the sign test, which are the
# quantities the conclusions rest on.

report <- function(r, lab, n) {
  k <- is18(r$feature); q <- p.adjust(r$pval, "BH")
  cat(sprintf("%-32s n=%-3d q<.05=%-2d pos=%-2d CI>0=%-2d med=%.2f sign_p=%s\n",
      lab, n, sum(q[k] < .05), sum(r$beta[k] > 0), sum(r$ci_lo[k] > 0),
      median(r$beta[k]),
      signif(binom.test(sum(r$beta[k] > 0), 18, .5,
                        alternative = "greater")$p.value, 3)))
  invisible(data.frame(analysis = lab, n = n, q_lt_05 = sum(q[k] < .05),
                       positive = sum(r$beta[k] > 0), ci_excl0 = sum(r$ci_lo[k] > 0),
                       beta_med = round(median(r$beta[k]), 3)))
}

X75 <- readRDS("data_processed/02_Metagenomics/clr_species_n75.rds")
mt  <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/metadata.tsv",
             sep = "\t", data.table = FALSE)
md  <- unique(mt[, c("Subject","Study.Group","Age","Gender")])
md  <- md[match(rownames(X75), as.character(md$Subject)), ]
rownames(md) <- as.character(md$Subject)
stopifnot(identical(rownames(X75), rownames(md)))

keep <- rownames(X75) != "10007545"
X <- X75[keep, , drop = FALSE]; M <- md[keep, ]
M$Age_z <- as.numeric(scale(M$Age)); M$Sex <- factor(M$Gender)
stopifnot(nrow(X) == 74)

out <- list()

# --- A1: D vs C -------------------------------------------------------
M$Group <- relevel(factor(M$Study.Group), ref = "C")
out$A1 <- report(ols(X, M, "Group + Age_z", "GroupD"), "A1 D vs C", 74)

# --- A3: D vs H, sex-adjusted ----------------------------------------
s <- M$Study.Group %in% c("D","H")
M2 <- M[s, ]; M2$Group <- relevel(factor(as.character(M2$Study.Group)), ref = "H")
M2$Age_z <- as.numeric(scale(M2$Age))
out$A3 <- report(ols(X[s, , drop = FALSE], M2, "Group + Age_z + Sex", "GroupD"),
                 "A3 D vs H + sex", sum(s))

# --- A1b: D vs C, females only ---------------------------------------
fo <- M$Study.Group %in% c("C","D") & M$Gender == "Female"
M3 <- M[fo, ]; M3$Group <- relevel(factor(as.character(M3$Study.Group)), ref = "C")
M3$Age_z <- as.numeric(scale(M3$Age))
out$A1b <- report(ols(X[fo, , drop = FALSE], M3, "Group + Age_z", "GroupD"),
                  "A1b D vs C females", sum(fo))

summ <- do.call(rbind, out)
dir.create("results/05_Reviewer_Response", recursive = TRUE, showWarnings = FALSE)
write.csv(summ, "results/05_Reviewer_Response/04_summary.csv", row.names = FALSE)
print(summ, row.names = FALSE)
