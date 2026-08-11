#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Isoform Infrastructure (Isopair Wrapper)
# Version: 7.0
#
# ⚠ REQUIRED — do not archive or delete. This is the ONLY producer of the
#   isoform-intrinsic infrastructure that 02_build_profiles_mashr.R consumes:
#   structures.rds, union_exons.rds, isoform_union_mapping.rds, cds.rds,
#   annotated_ue.rds and ptc.rds. `ptc.rds` in particular has no other producer
#   anywhere in the repository. It lived under archive/ until 2026-07-25, which
#   was wrong: nothing downstream can run without it.
#
#   Renamed 2026-07-25 from 02_build_profiles.R, which read as superseded by
#   02_build_profiles_mashr.R. It is not: the pair/profile outputs this script
#   also writes are NOT the published ones — those come from the mashr builder,
#   which additionally applies the 25% reference-share floor this script has no
#   concept of. The 01b_ prefix reflects the real order: it consumes
#   01_prepare_data_mashr.R's output and produces what 02_..._mashr.R needs.
#
# Purpose: Parse structures, build union exons, extract CDS, annotate region types
#          and compute PTC status. INFRASTRUCTURE ONLY — it builds no pairs and no
#          profiles. Those are 02_build_profiles_mashr.R's, and only its.
#
#          The infrastructure is isoform-intrinsic, so it is reused across changes in
#          sample scope. That reuse is safe ONLY while it still covers the expression
#          universe: the 2026-07-11 4-CT rescope regenerated expression_data.rds and
#          not this, leaving 3,323 of 95,623 isoforms (3.5%) with no structure and
#          therefore absent from the pair analysis, with no warning. That produced the
#          published 1,548 / 5,001 pair counts. 02 now refuses stale infrastructure
#          (its coverage guard); rerun THIS script whenever the universe changes.
#
# Input:
#   - RDS files from 01_prepare_data.R
#   - Isocall GTF (isoform structures)
#   - SQANTI corrected CDS GFF3 (CDS annotations for known + novel isoforms)
#
# Output (to output_dir/):
#   - structures.rds, union_exons.rds, isoform_union_mapping.rds
#   - cds.rds, annotated_ue.rds, ptc.rds
#
# Usage:
#   Rscript 01b_build_isoform_infrastructure.R [--data-dir DIR] [--output-dir DIR]

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
library(Isopair)
library(dplyr)

# ==============================================================================
# 0. Configuration
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_dir   <- file.path(.P$OUT, "data_mashr")
output_dir <- file.path(.P$OUT, "data_mashr")
if ("--data-dir" %in% args) {
  idx <- which(args == "--data-dir")
  if (idx < length(args)) data_dir <- args[idx + 1]
}
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  if (idx < length(args)) output_dir <- args[idx + 1]
}
# Rule in code: never write into a tree we read as input (R/load_config.R).
nmd_assert_writable(output_dir, .P, "--output-dir")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# GTF paths
isocall_gtf <- file.path(.P$DEPOSIT, "nmd_isocall_4ct.gtf.gz")
sqanti_cds_gff3 <- file.path(.P$SQANTI, "nmd_lungcells_corrected.cds.gff3")

cat("=== NMD Isoform Infrastructure (Isopair Wrapper) ===\n\n")

# ==============================================================================
# 1. Load prepared data
# ==============================================================================

# Only the expression matrix is needed, and only for its rownames: the isoform universe
# that structures are built over (section 2). gene_map, nmd_classification, dmso_samples
# and smg1i_samples were loaded here for the pair/profile loop removed 2026-07-26, and
# loading an input this script does not use is a false dependency -- it would fail a run
# that has infrastructure work to do but no classification yet.
#
# The 4-CT scope assertion that lived here went with them, deliberately and not by
# oversight: it checked nmd_classification, which this script no longer reads, and the
# stronger guard is in 02_build_profiles_mashr.R:79-89 where that data is actually used
# (26 samples, 13 DMSO / 13 Smg1i, no DD_ALI/DO_ALI/_DO_). Assert scope where it is
# consumed, not where it is merely in scope.
cat("Loading prepared data...\n")
expr_mat  <- readRDS(file.path(data_dir, "expression_data.rds"))

cat(sprintf("  Expression: %d isoforms x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# ==============================================================================
# 2. Isopair Infrastructure: structures, union exons, CDS, annotations
# ==============================================================================

# Parse isoform structures from isocall GTF
cat("\nParsing isoform structures...\n")
structures <- parseIsoformStructures(isocall_gtf,
  isoform_ids = rownames(expr_mat))
saveRDS(structures, file.path(output_dir, "structures.rds"))

# Build union exons
cat("Building union exons...\n")
ue <- buildUnionExons(structures)
saveRDS(ue$union_exons, file.path(output_dir, "union_exons.rds"))
saveRDS(ue$isoform_union_mapping, file.path(output_dir, "isoform_union_mapping.rds"))

# Extract CDS annotations from SQANTI corrected CDS GFF3
# This file has CDS features for both GENCODE and novel isoforms, using the
# same versioned IDs as isocall — no version stripping needed.
# Note: despite .gff3 extension, attributes are GTF format; specify format.
cat("Extracting CDS annotations...\n")
cat("  Loading SQANTI CDS GFF3...\n")
sqanti_gr <- rtracklayer::import(sqanti_cds_gff3, format = "gtf")
cds <- extractCdsAnnotations(sqanti_gr, isoform_ids = structures$isoform_id)

cat(sprintf("  CDS result: %d coding, %d unknown\n",
            sum(cds$coding_status == "coding"),
            sum(cds$coding_status == "unknown")))
saveRDS(cds, file.path(output_dir, "cds.rds"))

# Annotate region types
cat("Annotating region types...\n")
annotated_ue <- annotateRegionTypes(ue$isoform_union_mapping, cds)
saveRDS(annotated_ue, file.path(output_dir, "annotated_ue.rds"))

# Compute PTC status
cat("Computing PTC status...\n")
ptc <- computePtcStatus(structures, cds)
saveRDS(ptc, file.path(output_dir, "ptc.rds"))


# ==============================================================================
# 3. Done — infrastructure only, by design
# ==============================================================================
#
# The per-cell-type pair/profile loop that used to live here was REMOVED 2026-07-26.
# It wrote pairs_c{2,4}_{ct}.rds and profiles_c{2,4}_{ct}.rds -- the same filenames
# 02_build_profiles_mashr.R writes -- but WITHOUT the 25% reference-share floor, which
# exists only in 02. Two writers, one filename set, two different populations, and at
# least seven readers (05_final_report_mashr.Rmd, both gencode-scope reports,
# 06_orf_analysis_mashr.Rmd, 05r, 05k_b, predictor_comparison/01_extract_our_isoforms.R)
# that cannot tell which one they got. Which population reached Figures 3, 4 and 5 was
# decided by whichever script ran last, silently.
#
# The floored pairs are the published quantity -- the manuscript states the >= 25%
# criterion in the Figure 3 / Isopair sentence itself -- so this script's pairs were
# never the paper's, and in a correct run they were always overwritten. Removing them
# leaves exactly one writer, and makes this script's name true.
#
# Verified at removal: the six infrastructure outputs are byte-identical before and
# after (sha256), so nothing above this line changed.

cat("\n=== Infrastructure complete ===\n")
cat("Wrote to ", output_dir, ":\n", sep = "")
for (f in c("structures.rds", "union_exons.rds", "isoform_union_mapping.rds",
            "cds.rds", "annotated_ue.rds", "ptc.rds")) {
  fp <- file.path(output_dir, f)
  cat(sprintf("  %-28s %s\n", f,
              if (file.exists(fp)) sprintf("%.1f MB", file.size(fp) / 1024^2) else "*** MISSING ***"))
}
cat("\nNext: 02_build_profiles_mashr.R builds the pairs and profiles, applying the\n")
cat("      25% reference-share floor. It is the only writer of pairs_*/profiles_*.\n")
