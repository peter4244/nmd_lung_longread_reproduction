#!/usr/bin/env Rscript
# ==============================================================================
# 05k_b_utr5_refaug.R
#
# Ref-AUG-projected 5'UTR feature scan — sibling to 05k_utr5_all_isoforms.R.
#
# Motivation: 05k_utr5_all_isoforms.R uses the comparator's TD2-called CDS to
# define the 5'UTR boundary. TD2 has documented PTC-avoidance bias (Figure 3),
# so for PTC-containing comparators TD2's chosen start sits DOWNSTREAM of the
# true biological start. That inflates the apparent 5'UTR and demotes true
# protein-coding sequence to "uORF". This script re-scans the 5'UTR using the
# ref-AUG-projected position in the comparator (from ref_atg_analysis.rds) as
# the boundary instead.
#
# Scope: all comparator isoforms appearing in ref_atg_analysis.rds$c2 or $c4
# with category ∈ {effectively_ptc, no_downstream_ejc, truncated_no_ejc}
# (categories where a ref-AUG-projected position exists). This is broader than
# the §4 ENST-reference scope so the output is reusable for other analyses;
# downstream consumers filter to their own scope.
#
# Inputs:
#   - data_mashr/structures.rds
#   - data_mashr/analysis_cache/ref_atg_analysis.rds
#   - data_mashr/profiles_c2_allsamples.rds (for set/group labels)
#   - data_mashr/profiles_c4_allsamples.rds
#   - FASTA: sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta
#
# Output:
#   - data_mashr/analysis_cache/utr5_features_refaug.rds
#
# Output schema matches utr5_features_all.rds$isoform_features so downstream
# consumers can swap inputs. Adds a `boundary_source = "ref_aug_projected"`
# metadata field. Original utr5_features_all.rds (TD2-anchored) is NOT modified.
#
# Pre-registered in figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md §4.1.
# Usage:
#   Rscript 05k_b_utr5_refaug.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(Isopair)
})

# setwd(.P$ISOPAIR) REMOVED 2026-07-28 -- same defect as 05k_utr5_all_isoforms.R:30.
# .P$ISOPAIR is the legacy authoring tree and does not exist for a reader; setwd() to a
# missing directory is a hard error. Measured: this file has 0 relative read/source calls,
# so the call was dead weight and the only thing stopping the step running. The graph could
# not see it -- setwd() is neither a read nor a write, so needs_update scored FALSE. W80.

# --- DATA_DIR: the half-repoint this script used to have ------------------------------
# Every data read below was a bare relative path under the setwd() above, i.e. the LEGACY
# tree, while every write already went to .P$OUT. That combination silently builds the
# analysis cache from the PUBLISHED pairs and files it beside a rebuild, so a downstream
# report mixes vintages and looks fine. Anchor the reads too. (2026-07-26)
#
# Default is .P$OUT/data_mashr. Note that 01b and 02 write the SAME filenames there and 02
# overwrites 01b -- only 02's are floored, and only the floored ones are the published
# quantity. Point --data-dir at the floored tree explicitly rather than trusting the default.
DATA_DIR <- file.path(.P$OUT, "data_mashr")
local({
  a <- commandArgs(trailingOnly = TRUE)
  if ("--data-dir" %in% a) { i <- which(a == "--data-dir"); if (i < length(a)) DATA_DIR <<- a[i + 1] }
})
DATA_DIR <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(DATA_DIR)) {
  stop(sprintf("DATA_DIR does not exist: %s\n  Run 01_prepare_data_mashr.R -> 01b -> 02 first, or pass --data-dir.", DATA_DIR))
}
cat("data dir: ", DATA_DIR, "\n", sep = "")
# Rule in code: this script WRITES here, so it must not be an input root. A --data-dir
# pointing at a reference tree would overwrite the copies verification compares against.
nmd_assert_writable(DATA_DIR, .P, "--data-dir")


FASTA_PATH  <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")
OUTPUT_PATH <- file.path(DATA_DIR, "analysis_cache", "utr5_features_refaug.rds")   # the cache belongs beside the data it was built FROM   # anchored: setwd() is into an INPUT root
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

cat("=== Ref-AUG-Projected 5'UTR Feature Scan ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ─────────────────────────────────────────────────────────────────────
# 1. Load ref_atg analysis + pair-level metadata for scope and labels
# ─────────────────────────────────────────────────────────────────────
ref_atg <- readRDS(file.path(DATA_DIR, "analysis_cache/ref_atg_analysis.rds"))
structures <- readRDS(file.path(DATA_DIR, "structures.rds"))
profiles_c2 <- readRDS(file.path(DATA_DIR, "profiles_c2_allsamples.rds"))
profiles_c4 <- readRDS(file.path(DATA_DIR, "profiles_c4_allsamples.rds"))

TRACEABLE <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")

# Filter to ref-AUG-traceable comparator pairs
ra_c2 <- ref_atg$c2[ref_atg$c2$category %in% TRACEABLE, ]
ra_c4 <- ref_atg$c4[ref_atg$c4$category %in% TRACEABLE, ]
cat(sprintf("Ref-AUG-traceable pairs in cache: c2 = %d (NMD), c4 = %d (Control)\n",
            nrow(ra_c2), nrow(ra_c4)))

# Unique comparator isoforms to scan (a comparator may appear in multiple
# pairs only via different references; we scan each comparator once)
target_ids <- unique(c(ra_c2$comparator_isoform_id, ra_c4$comparator_isoform_id))
cat(sprintf("Unique comparator isoforms in scope: %d\n", length(target_ids)))

# Retain ref_atg_genomic and comp_stop_tx_pos per isoform.
# Prefer the c2 row (NMD-side) if an isoform appears on both sides; we expect
# this not to happen for comparators (an isoform is NMD or non-NMD, not both),
# but defensive dedup keeps the logic safe.
extract_traceable <- function(df, set_label) {
  data.frame(
    isoform_id      = df$comparator_isoform_id,
    ref_atg_genomic = df$ref_atg_genomic,
    comp_stop_tx_pos = df$comp_stop_tx_pos,
    strand          = df$strand,
    category        = df$category,
    set             = set_label,
    stringsAsFactors = FALSE
  )
}
scope_dt <- rbind(extract_traceable(ra_c2, "NMD"),
                  extract_traceable(ra_c4, "Control"))
scope_dt <- scope_dt[!duplicated(scope_dt$isoform_id), ]
cat(sprintf("Per-isoform scope (after dedup): %d isoforms\n", nrow(scope_dt)))
cat("  Set composition:\n"); print(table(scope_dt$set))
cat("  Category composition:\n"); print(table(scope_dt$category))

# ─────────────────────────────────────────────────────────────────────
# 2. Compute ref-AUG-traced stop genomic position per isoform.
# `comp_stop_tx_pos` lives in transcript space; convert to genomic via
# Isopair::transcriptToGenomic.
# ─────────────────────────────────────────────────────────────────────
struct_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i)
    list(starts = structures$exon_starts[[i]],
         ends   = structures$exon_ends[[i]])),
  structures$isoform_id)

cat(sprintf("\nConverting comp_stop_tx_pos to genomic coords for %d isoforms...\n",
            nrow(scope_dt)))
stop_genomic <- rep(NA_integer_, nrow(scope_dt))
for (i in seq_len(nrow(scope_dt))) {
  iso <- scope_dt$isoform_id[i]
  cs  <- struct_lookup[[iso]]
  if (is.null(cs)) next
  tp  <- scope_dt$comp_stop_tx_pos[i]
  if (is.na(tp)) next
  stop_genomic[i] <- Isopair::transcriptToGenomic(
    tp, cs$starts, cs$ends, scope_dt$strand[i])
}
scope_dt$stop_genomic <- stop_genomic
cat(sprintf("  Successfully mapped to genomic: %d / %d\n",
            sum(!is.na(stop_genomic)), nrow(scope_dt)))

# ─────────────────────────────────────────────────────────────────────
# 3. Build synthetic cds_metadata with ref-AUG-projected coordinates.
# cds.rds convention: cds_start = smaller genomic coord, cds_stop = larger,
# strand carries reading direction. scan5UtrFeatures expects this.
# ─────────────────────────────────────────────────────────────────────
scope_dt$cds_start <- pmin(scope_dt$ref_atg_genomic, scope_dt$stop_genomic)
scope_dt$cds_stop  <- pmax(scope_dt$ref_atg_genomic, scope_dt$stop_genomic)
ok <- !is.na(scope_dt$cds_start) & !is.na(scope_dt$cds_stop)
scope_dt <- scope_dt[ok, ]
cat(sprintf("Isoforms with both ref-AUG and traced stop genomic positions: %d\n",
            nrow(scope_dt)))

cds_meta_refaug <- data.frame(
  isoform_id    = scope_dt$isoform_id,
  coding_status = "coding",
  cds_start     = as.integer(scope_dt$cds_start),
  cds_stop      = as.integer(scope_dt$cds_stop),
  orf_length    = NA_integer_,
  strand        = scope_dt$strand,
  stringsAsFactors = FALSE
)

# Structures filtered to in-scope isoforms
str_scope <- structures[structures$isoform_id %in% scope_dt$isoform_id, ]

# ─────────────────────────────────────────────────────────────────────
# 4. Extract sequences from FASTA
# ─────────────────────────────────────────────────────────────────────
id_file <- tempfile(fileext = ".txt")
writeLines(scope_dt$isoform_id, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH)
t_awk <- proc.time()
raw_output <- system(awk_cmd, intern = TRUE)
cat(sprintf("\nawk completed in %.1f seconds\n", (proc.time() - t_awk)[3]))
unlink(id_file)
split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
seq_vec <- setNames(
  toupper(vapply(split_lines, `[`, character(1), 2)),
  vapply(split_lines, `[`, character(1), 1))
n_with_seq <- sum(scope_dt$isoform_id %in% names(seq_vec))
cat(sprintf("Sequences loaded: %d  (in-scope: %d / %d)\n",
            length(seq_vec), n_with_seq, nrow(scope_dt)))

# ASSERT WHAT THE EXTRACTION DID, NOT THAT THE CALL RETURNED (W122). system(intern = TRUE)
# raises a WARNING and never an error on non-zero exit, so a failed awk returns character(0),
# strsplit/vapply reduce that to an EMPTY seq_vec without complaint, and scan5UtrFeatures then
# runs over no sequences and writes a plausible-looking empty feature table at exit 0. A
# PARTIAL shortfall is legitimate and is reported above; zero is not.
.awk_status <- attr(raw_output, "status")
if (!is.null(.awk_status) && .awk_status != 0)
  stop(sprintf("awk sequence extraction exited %s; no sequences were extracted from %s",
               .awk_status, FASTA_PATH))
if (nrow(scope_dt) > 0 && n_with_seq == 0)
  stop(sprintf("awk returned no sequence for ANY of the %d in-scope isoforms in %s",
               nrow(scope_dt), FASTA_PATH))

# ─────────────────────────────────────────────────────────────────────
# 5. Scan 5'UTR features using ref-AUG-projected boundaries
# ─────────────────────────────────────────────────────────────────────
cat("\nScanning ref-AUG-anchored 5'UTR features...\n")
t_scan <- proc.time()
isoform_features <- scan5UtrFeatures(
  structures   = str_scope,
  cds_metadata = cds_meta_refaug,
  sequences    = seq_vec,
  min_orf_nt   = 9L,
  verbose      = TRUE)
elapsed <- round((proc.time() - t_scan)[3], 1)
cat(sprintf("Completed in %.1f seconds\n", elapsed))
cat(sprintf("Features computed for %d isoforms\n", nrow(isoform_features)))

# ATG validation
n_validated <- sum(!is.na(isoform_features$atg_validated))
n_pass <- sum(isoform_features$atg_validated, na.rm = TRUE)
cat(sprintf("ATG validation: %d / %d (%.1f%%) pass\n",
            n_pass, n_validated,
            if (n_validated > 0) 100 * n_pass / n_validated else NA))

# Exclusion flag (matches 05k_utr5_all_isoforms.R logic)
isoform_features$excluded <- isoform_features$utr5_length == 0L |
  (!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)
n_excluded <- sum(isoform_features$excluded)
n_included <- sum(!isoform_features$excluded)
cat(sprintf("Excluded: %d  Included: %d\n", n_excluded, n_included))

# Carry over the set/category labels
scope_lookup <- setNames(scope_dt$set,      scope_dt$isoform_id)
cat_lookup   <- setNames(scope_dt$category, scope_dt$isoform_id)
isoform_features$set      <- scope_lookup[isoform_features$isoform_id]
isoform_features$category <- cat_lookup[isoform_features$isoform_id]

cat("\nSet x category counts:\n")
print(table(isoform_features$set, isoform_features$category, useNA = "ifany"))

# ─────────────────────────────────────────────────────────────────────
# 6. Save
# ─────────────────────────────────────────────────────────────────────
output <- list(
  isoform_features = isoform_features,
  metadata = list(
    run_timestamp        = Sys.time(),
    script_path          = "05k_b_utr5_refaug.R",
    fasta_path           = FASTA_PATH,
    boundary_source      = "ref_aug_projected",
    scope_description    = "ref_atg_analysis.rds c2 ∪ c4 ∩ traceable categories",
    n_in_scope           = nrow(scope_dt),
    n_with_sequence      = n_with_seq,
    n_features_computed  = nrow(isoform_features),
    n_excluded           = n_excluded,
    n_included           = n_included,
    elapsed_scan_seconds = elapsed,
    sibling_note         = "Companion to utr5_features_all.rds (TD2-anchored). Same schema."
  )
)

saveRDS(output, OUTPUT_PATH)
cat(sprintf("\nSaved to: %s\n", OUTPUT_PATH))
cat("Finished:", format(Sys.time()), "\n")
