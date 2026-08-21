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
  # every coefficient the model returns, so BH below spans the full family
  terms <- rownames(B)[-1]
  all <- do.call(rbind, lapply(terms, function(k) {
    se <- sqrt(iv[k, k] * s2); b <- B[k, ]
    data.frame(feature = colnames(X), term = k, beta = b, se = se,
               ci_lo = b - qt(.975, dfr) * se, ci_hi = b + qt(.975, dfr) * se,
               pval = 2 * pt(-abs(b / se), dfr), df = dfr, row.names = NULL)
  }))
  all$qval <- p.adjust(all$pval, "BH")
  attr(all, "family_size") <- nrow(all)
  out <- all[all$term == term, ]
  attr(out, "family_size") <- nrow(all)
  out
}

# Multiple-testing scope
# -----------------------
# BH is applied across every coefficient the model returns, matching the
# family used by MaAsLin2 in the primary model and verified empirically in
# 03_verify_BH_scope.R. Each model here returns three coefficients per
# feature, so the family is 3 x 4,795 = 14,385 tests, the same denominator
# as the primary model. Correcting within a single contrast (4,795) would
# use a narrower family than the analysis it is compared against.

report <- function(r, lab, n) {
  k <- is18(r$feature); q <- r$qval
  # No sign test: with 18 of 18 positive the p-value is 0.5^18 by construction
  # and carries no information about the data. Direction counts and intervals
  # are reported instead.
  cat(sprintf("%-32s n=%-3d family=%-6d q<.05=%-2d pos=%-2d CI>0=%-2d med=%.2f\n",
      lab, n, attr(r, "family_size"), sum(q[k] < .05),
      sum(r$beta[k] > 0), sum(r$ci_lo[k] > 0), median(r$beta[k])))
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
