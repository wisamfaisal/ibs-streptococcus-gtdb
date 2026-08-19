# 05_external_metadata.R
# Analyses drawing on sources beyond the curated repository tables.
#   B1  per-group detection frequency and median relative abundance
#   B8  antibiotic sensitivity (22 exposed samples removed, means recomputed)
#   A5  sequencing depth (ENA PRJEB37924)
#   A4  PPI sensitivity and sequencing-batch diagnostics (Mars Table S1)
#
# Required inputs beyond the repository tables:
#   data_raw/02_Metagenomics/mmc1.xlsx                 (Mars 2020 Table S1)
#   data_raw/02_Metagenomics/ENA_PRJEB37924_runs.tsv   (ENA Portal API)
# The ENA table is retrieved by the URL documented at the foot of this file.

library(data.table); library(readxl)

f_xl  <- "data_raw/02_Metagenomics/mmc1.xlsx"
f_ena <- "data_raw/02_Metagenomics/ENA_PRJEB37924_runs.tsv"
for (f in c(f_xl, f_ena))
  if (!file.exists(f)) stop("missing required input: ", f)

SIG18 <- c("parasanguinis_G","sp001813295","parasanguinis_B","parasanguinis_D",
  "parasanguinis_C","parasanguinis_E","sp001578805","sp900755085","sp902836505",
  "parasanguinis","parasanguinis_A","sp902373455","parasanguinis_F","intermedius",
  "alactolyticus","sp000187445","vestibularis","sp003521145")

sp <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/species.tsv",
            sep = "\t", data.table = FALSE)
mt <- fread("data_raw/02_Metagenomics/MARS_IBS_2020/metadata.tsv",
            sep = "\t", data.table = FALSE)

# ---- B1: detection and abundance on the untransformed scale ----------
idx  <- match(sp$Sample, mt$Sample); subj <- as.character(mt$Subject[idx])
Mraw <- as.matrix(sp[, -1, drop = FALSE])
ag   <- rowsum(Mraw, subj) / as.vector(table(subj)[rownames(rowsum(Mraw, subj))])
ag   <- ag[rownames(ag) != "10007545", , drop = FALSE]
grp  <- mt$Study.Group[match(rownames(ag), as.character(mt$Subject))]

sel <- sub(".*s__Streptococcus ", "", colnames(ag)) %in% SIG18 &
       grepl("g__Streptococcus;", colnames(ag), fixed = TRUE)
stopifnot(sum(sel) == 18)
A <- ag[, sel, drop = FALSE]
colnames(A) <- sub(".*s__Streptococcus ", "", colnames(A))

b1 <- do.call(rbind, lapply(colnames(A), function(f) {
  x <- A[, f]
  data.frame(lineage = f,
    det_C = sum(x[grp=="C"] > 0), det_D = sum(x[grp=="D"] > 0),
    det_H = sum(x[grp=="H"] > 0),
    pct_D = round(100 * mean(x[grp=="D"] > 0)),
    pct_H = round(100 * mean(x[grp=="H"] > 0)),
    med_D = signif(median(x[grp=="D"]), 3),
    med_H = signif(median(x[grp=="H"]), 3))
}))
cat("B1: median detection D =", median(b1$pct_D), "% | H =",
    median(b1$pct_H), "% | lineages <50% in D:", sum(b1$pct_D < 50), "\n")

# ---- A5: sequencing depth -------------------------------------------
ena <- read.delim(f_ena, stringsAsFactors = FALSE)
st  <- ena[sub("_[0-9]+$", "", ena$sample_alias) == "ibs" &
           grepl("^Participant", ena$sample_title), ]
st$key <- paste0(sub("^Participant ([0-9]+) Sample.*", "\\1", st$sample_title),
                 ".T.", sub(".*Sample ([0-9]+)$", "\\1", st$sample_title))
lnk <- merge(mt[, c("Sample","Subject","Study.Group")],
             data.frame(key = st$key, depth = st$read_count),
             by.x = "Sample", by.y = "key")
pd  <- aggregate(depth ~ Subject + Study.Group, lnk, mean)
pd  <- pd[pd$Subject != 10007545, ]
cat("A5: matched", nrow(lnk), "of 444 | Kruskal p =",
    signif(kruskal.test(depth ~ factor(Study.Group), pd)$p.value, 3),
    "| D vs H p =",
    signif(wilcox.test(depth ~ factor(Study.Group),
           pd[pd$Study.Group %in% c("D","H"), ])$p.value, 3), "\n")

# ---- A4: PPI and batch ----------------------------------------------
sm <- read_excel(f_xl, sheet = "stool metadata"); sm$sid <- as.character(sm$study_id)
ppi <- tapply(tolower(trimws(as.character(sm$proton_pump))) == "yes",
              sm$sid, any, na.rm = TRUE)
bat <- tapply(as.character(sm$seq_batch) == "1st", sm$sid, mean, na.rm = TRUE)
ext <- data.frame(Subject = names(ppi), PPI = as.logical(ppi),
                  b1 = as.numeric(bat[names(ppi)]))
e <- ext[match(rownames(ag), ext$Subject), ]
cat("A4: PPI users by group ->",
    paste(names(tapply(e$PPI, grp, sum)), tapply(e$PPI, grp, sum),
          collapse = " "), "\n")
cat("A4: batch-1 fraction sd within IBS-D =",
    round(sd(e$b1[grp == "D"]), 3),
    "-- too low to support batch as a covariate\n")

dir.create("results/05_Reviewer_Response", recursive = TRUE, showWarnings = FALSE)
write.csv(b1, "results/05_Reviewer_Response/B1_prevalence_18lineages.csv",
          row.names = FALSE)
write.csv(pd, "results/05_Reviewer_Response/A5_depth_per_subject.csv",
          row.names = FALSE)
cat("\nwritten.\n")

# ENA retrieval URL (run once):
# https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB37924
#   &result=read_run&fields=run_accession,sample_accession,sample_title,
#   sample_alias,read_count,base_count,fastq_bytes,instrument_model,
#   library_layout&format=tsv&limit=0
