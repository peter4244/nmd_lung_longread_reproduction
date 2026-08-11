#!/usr/bin/env Rscript
# =============================================================================
# 00b_build_cds_exons.R — producer for data_mashr/cds_exons.rds
#
# WHY THIS EXISTS
# ---------------
# `cds_exons.rds` is read by three kept pipeline files:
#   - 02_build_profiles_mashr.R:99   (copies it into data_mashr/)
#   - 03b_rebuild_cache.R:54         (-> compareIsoformFrames(profiles, cds_exons, cds))
#   - 05_final_report_mashr.Rmd:383
# but until 2026-07-21 NOTHING in the repo produced it. The on-disk copy is dated
# Mar 12 12:28, a day after its six sibling infrastructure objects (Mar 11 22:59),
# i.e. it was made by an ad-hoc interactive step that was never scripted. Under the
# project's reproducibility standard (raw reads + code, no interim files deposited)
# that was an unreproducible link in the §4 chain.
#
# VERIFIED: the output of this script is IDENTICAL to the existing
# data_mashr/cds_exons.rds — 755,368 rows x 5 cols, same column names, and
# all.equal() TRUE on content after sorting by (isoform_id, cds_exon_rank).
#
# Run AFTER the infrastructure build (structures.rds must exist) and BEFORE
# 02_build_profiles_mashr.R.
# =============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(Isopair)
  library(rtracklayer)
})

# --- DATA_DIR and the CDS source, both repointed 2026-07-26 --------------------------
# This read structures.rds from the LEGACY data_mashr and wrote cds_exons.rds straight back
# into it -- a write into an input root, which tools/check_output_roots.R exists to prevent.
# It also took the CDS GFF3 from the legacy sqanti tree, though the identical file is
# deposited and hash-verified; 01b_build_isoform_infrastructure.R already reads the
# deposited copy. Both now resolve through config.
#
# Pass --data-dir explicitly: 01b and 02 write the same filenames to the same directory and
# only 02's pairs carry the 25% reference-share floor.
DATA_DIR <- file.path(.P$OUT, "data_mashr")
local({
  a <- commandArgs(trailingOnly = TRUE)
  if ("--data-dir" %in% a) { i <- which(a == "--data-dir"); if (i < length(a)) DATA_DIR <<- a[i + 1] }
})
DATA_DIR <- normalizePath(DATA_DIR, mustWork = FALSE)

# Create the parent before writing. This step SURVIVED the 2026-07-28 clean-room run only
# because 01_prepare_data_mashr.R happened to create data_mashr/ at index 20, twelve steps
# earlier -- an accident of run_order, not a guarantee. Its sibling 05s_orfik_scan.R had no
# such luck and lost 75 minutes of completed work to the same missing directory. W81.
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DATA_DIR)) {
  stop(sprintf("DATA_DIR does not exist: %s\n  Run 01_prepare_data_mashr.R -> 01b first, or pass --data-dir.", DATA_DIR))
}
cat("data dir: ", DATA_DIR, "\n", sep = "")
# Rule in code: this script WRITES here, so it must not be an input root. A --data-dir
# pointing at a reference tree would overwrite the copies verification compares against.
nmd_assert_writable(DATA_DIR, .P, "--data-dir")

# Same CDS source the rest of the infrastructure uses, now the DEPOSITED copy (in
# MANIFEST.sha256). Despite the .gff3 extension the attributes are GTF-format, hence
# format = "gtf".
SQANTI_CDS_GFF3 <- file.path(.P$SQANTI, "nmd_lungcells_corrected.cds.gff3")

stopifnot(file.exists(SQANTI_CDS_GFF3))
structures <- readRDS(file.path(DATA_DIR, "structures.rds"))

cat("Importing SQANTI CDS annotations...\n")
sqanti_gr <- rtracklayer::import(SQANTI_CDS_GFF3, format = "gtf")

cat("Extracting per-exon CDS segments...\n")
cds_exons <- extractCdsExons(sqanti_gr,
                             isoform_ids = structures$isoform_id,
                             verbose     = TRUE)

cat(sprintf("cds_exons: %d rows x %d cols\n", nrow(cds_exons), ncol(cds_exons)))
saveRDS(cds_exons, file.path(DATA_DIR, "cds_exons.rds"))
cat("Wrote", file.path(DATA_DIR, "cds_exons.rds"), "\n")
