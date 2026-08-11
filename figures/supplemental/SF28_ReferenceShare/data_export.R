#!/usr/bin/env Rscript
# SF28 — Reference-isoform share of gene expression (post 25% reference-share floor).
#
# Rebuild of the reference-share export (originally in the deleted
# PairSetDescriptives/data_export.R; recover with
#   git show f6d96bc^:figures/SupplementalFigures/PairSetDescriptives/data_export.R
# ). Two changes vs the original for consistency with the 2026-07-10 floor:
#   1. Denominator = ALL isoforms of the gene (not just non-NMD), matching the
#      floor's all-isoform basis.
#   2. DMSO basis = dmso_samples[["all_samples"]] (the basis the reference
#      selection + floor use in 02_build_profiles_mashr.R), not the 4-CT union.
# Population = the FLOORED pop_BC (profiles are already floored), so every gene
# here has ref_fraction_of_gene >= 0.25 by construction (retained-only view).
#
# Output (SF28_ReferenceShare/data/):
#   ref_expression_fraction.tsv  (gene_id, reference_isoform_id, ref_fraction_of_gene)
#   descriptives_summary.tsv     (metric, value)

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressMessages(library(data.table))

# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this

# exported panels from the 2026-03-11 structure universe while the reports use the

# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.

DM  <- nmd_data_dir(.P)
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)))),
          error = function(e) ".")
OUT <- file.path(HERE, "data"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

expr        <- readRDS(file.path(DM, "expression_data.rds"))
structures  <- as.data.table(readRDS(file.path(DM, "structures.rds")))[, .(isoform_id, gene_id)]
dmso        <- readRDS(file.path(DM, "dmso_samples.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

# pop_BC = gene+reference present in both C2 and C4 (already floored in 02).
shared <- merge(unique(profiles_c2[, .(gene_id, reference_isoform_id)]),
                unique(profiles_c4[, .(gene_id, reference_isoform_id)]),
                by = c("gene_id", "reference_isoform_id"))

# all-isoform gene totals on the all_samples DMSO basis (floor's basis).
dmso_all   <- dmso[["all_samples"]]
stopifnot(all(dmso_all %in% colnames(expr)))
dmso_means <- rowMeans(expr[, dmso_all, drop = FALSE])
gm <- structures[isoform_id %in% names(dmso_means)]
gm[, mean_dmso := dmso_means[isoform_id]][is.na(mean_dmso), mean_dmso := 0]
gene_total <- gm[, .(gene_total = sum(mean_dmso)), by = gene_id]

fra <- merge(shared, gene_total, by = "gene_id", all.x = TRUE)
fra[, ref_dmso := dmso_means[reference_isoform_id]][is.na(ref_dmso), ref_dmso := 0]
fra[, ref_fraction_of_gene := ifelse(gene_total > 0, ref_dmso / gene_total, NA_real_)]
fra <- fra[!is.na(ref_fraction_of_gene)]

# Sanity: floored population — every retained gene is >= 25% by construction.
stopifnot(min(fra$ref_fraction_of_gene) >= 0.25 - 1e-9)

fwrite(fra[, .(gene_id, reference_isoform_id, ref_fraction_of_gene)],
       file.path(OUT, "ref_expression_fraction.tsv"), sep = "\t")

pct <- fra$ref_fraction_of_gene * 100
summary_dt <- data.table(
  metric = c("n_genes_in_pop_BC", "median_ref_share_pct", "pct_ge_50", "floor_pct"),
  value  = c(nrow(fra), round(median(pct), 1), round(mean(pct >= 50) * 100, 1), 25))
fwrite(summary_dt, file.path(OUT, "descriptives_summary.tsv"), sep = "\t")

cat(sprintf("SF28: n=%d genes  median=%.1f%%  >=50%%: %.1f%%  (floor 25%%, all-iso, all_samples DMSO)\n",
            nrow(fra), median(pct), mean(pct >= 50) * 100))
