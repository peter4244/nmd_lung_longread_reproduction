#!/usr/bin/env Rscript
# SF27 — expressed non-NMD isoform count per pop_BC gene (post 25% floor).
#
# Standalone export (no prior producer on disk). Reproduces the canonical Rmd
# §1 `sec1-iso-per-gene` definition: for each pop_BC gene, count isoforms that
# are (a) in structures, (b) expressed (in expression_data), (c) NOT classified
# NMD-susceptible in the all_samples analysis. Population = the FLOORED pop_BC
# (profiles already floored). See REFERENCE_FLOOR_SF_INVENTORY.md.
#
# Output (SF27_IsoformsPerGene/data/): isoforms_per_gene.tsv (gene_id, n_isoforms),
#   descriptives_summary.tsv (metric, value).

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

expr        <- readRDS(file.path(DM, "expression_data.rds"))
structures  <- as.data.table(readRDS(file.path(DM, "structures.rds")))[, .(isoform_id, gene_id)]
nmd_class   <- readRDS(file.path(DM, "nmd_classification.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

mk <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
popbc <- unique(profiles_c2[mk(profiles_c2) %in% intersect(mk(profiles_c2), mk(profiles_c4)), gene_id])

nmd_set <- nmd_class[["all_samples"]]$nmd
sec1 <- structures[gene_id %in% popbc &
                   isoform_id %in% rownames(expr) &
                   !isoform_id %in% nmd_set]
ipg <- sec1[, .(n_isoforms = .N), by = gene_id]
# genes with 0 non-NMD expressed isoforms keep a 0 row (rare; preserves denominator)
zero <- setdiff(popbc, ipg$gene_id)
if (length(zero)) ipg <- rbind(ipg, data.table(gene_id = zero, n_isoforms = 0L))

fwrite(ipg[order(gene_id)], file.path(OUT, "isoforms_per_gene.tsv"), sep = "\t")
summ <- data.table(
  metric = c("n_genes_in_pop_BC", "median_isoforms_per_gene", "iqr_isoforms_per_gene"),
  value  = c(nrow(ipg), as.integer(median(ipg$n_isoforms)), as.integer(IQR(ipg$n_isoforms))))
fwrite(summ, file.path(OUT, "descriptives_summary.tsv"), sep = "\t")
cat(sprintf("SF27: %d genes  median=%d isoforms/gene  IQR=%d\n",
            nrow(ipg), as.integer(median(ipg$n_isoforms)), as.integer(IQR(ipg$n_isoforms))))
