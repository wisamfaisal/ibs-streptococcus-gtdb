# 03_verify_BH_scope.R
# Establishes empirically that MaAsLin2 applied BH across all 14,385
# returned coefficients (2 group contrasts + age, x 4,795 features),
# not across the 9,590 group-contrast tests alone.

pub <- read.delim("results/03_Metagenomics_FINAL_model/all_results.tsv",
                  stringsAsFactors = FALSE)
stopifnot(nrow(pub) == 14385, all(table(pub$value) == 4795))

# candidate scope A: all returned coefficients
qA <- p.adjust(pub$pval, "BH")
dA <- max(abs(qA - pub$qval))

# candidate scope B: group contrasts only
g  <- pub$value != "Age"
qB <- p.adjust(pub$pval[g], "BH")
dB <- max(abs(qB - pub$qval[g]))

cat("scope = all 14,385   -> max |diff| =", format(dA, scientific = TRUE), "\n")
cat("scope = 9,590 only   -> max |diff| =", format(dB, scientific = TRUE), "\n")
stopifnot(dA < 1e-12, dB > 1e-3)

sig <- pub[pub$qval < 0.05, ]
stopifnot(nrow(sig) == 50, sum(sig$value != "Age") == 49)
cat("passing coefficients:", nrow(sig),
    "| group-level:", sum(sig$value != "Age"),
    "| age term:", sum(sig$value == "Age"), "\n")
cat("VERIFIED: BH family = 14,385\n")
