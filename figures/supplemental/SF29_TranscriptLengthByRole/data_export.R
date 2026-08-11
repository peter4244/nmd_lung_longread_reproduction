#!/usr/bin/env Rscript
# SF29 — transcript length by role in the pop_BC triplet (post 25% floor).
#
# Standalone export reproducing the canonical Rmd §1 `sec1-tx-length` definition:
# transcript length = sum(exon lengths) per isoform; one gene-matched triplet
# (Reference, NMD comparator, Control comparator) per pop_BC gene; paired
# Wilcoxon across roles. Population = FLOORED pop_BC. See REFERENCE_FLOOR_SF_INVENTORY.md.
#
# Output (SF29_TranscriptLengthByRole/data/):
#   tx_length_by_role_long.tsv (role, length_nt), descriptives_summary.tsv (metric, value).

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressMessages(library(data.table))
# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this
# exported panels from the 2026-03-11 structure universe while the reports use the
# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.
DM <- nmd_data_dir(.P)
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)))),
          error = function(e) ".")
OUT <- file.path(HERE, "data"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

structures  <- readRDS(file.path(DM, "structures.rds"))
p2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
p4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

tx_len <- setNames(
  vapply(seq_len(nrow(structures)),
         function(i) sum(structures$exon_ends[[i]] - structures$exon_starts[[i]] + 1L),
         integer(1)),
  structures$isoform_id)

trip <- merge(
  p2[, .(gene_id, reference_isoform_id, comparator_isoform_id_nmd = comparator_isoform_id)],
  p4[, .(gene_id, reference_isoform_id, comparator_isoform_id_ctrl = comparator_isoform_id)],
  by = c("gene_id", "reference_isoform_id"))
trip[, ref_len  := tx_len[reference_isoform_id]]
trip[, nmd_len  := tx_len[comparator_isoform_id_nmd]]
trip[, ctrl_len := tx_len[comparator_isoform_id_ctrl]]
trip <- trip[complete.cases(trip[, .(ref_len, nmd_len, ctrl_len)])]

long <- rbindlist(list(
  data.table(role = "Reference",          length_nt = trip$ref_len),
  data.table(role = "NMD comparator",     length_nt = trip$nmd_len),
  data.table(role = "Control comparator", length_nt = trip$ctrl_len)))
fwrite(long, file.path(OUT, "tx_length_by_role_long.tsv"), sep = "\t")

pw <- function(a, b) wilcox.test(a, b, paired = TRUE, exact = FALSE)$p.value
summ <- data.table(
  metric = c("n_pop_bc_triplets", "median_ref_length_nt", "median_nmd_comp_length_nt",
             "median_control_comp_length_nt", "paired_wilcox_p_nmd_vs_ref",
             "paired_wilcox_p_nmd_vs_ctrl", "paired_wilcox_p_ctrl_vs_ref"),
  value  = c(nrow(trip), median(trip$ref_len), median(trip$nmd_len), median(trip$ctrl_len),
             signif(pw(trip$nmd_len, trip$ref_len), 3),
             signif(pw(trip$nmd_len, trip$ctrl_len), 3),
             signif(pw(trip$ctrl_len, trip$ref_len), 3)))
fwrite(summ, file.path(OUT, "descriptives_summary.tsv"), sep = "\t")
cat(sprintf("SF29: %d triplets  medians ref/nmd/ctrl = %d/%d/%d nt\n",
            nrow(trip), as.integer(round(median(trip$ref_len))),
            as.integer(round(median(trip$nmd_len))), as.integer(round(median(trip$ctrl_len)))))
