#!/usr/bin/env Rscript
# Data export for new Fig 3 (Isopair + PTC combined multipanel).
#
# 2026-06-13 PM rewrite — population structure based on Pete's clarified scope:
#
# Scope policy:
#   "TD2 CDS annotations are unreliable, so for analyses that depend on
#    identifying the stop codon, we limit to isoform pairs where reference
#    AUG tracing can be performed."
#
# Three populations used in this figure:
#
#   pop_BC (size from the data; 1,548 published, 1,578 in the rebuild):
#       Stage 2 gene-matched pairs. No CDS dependency. Used by Panels B
#       (sequence similarity) and C (splice event prevalence).
#
#   pop_traceable (n=2,289 NMD / 1,763 Control):
#       Pairs where ref-AUG tracing produced a valid ORF position
#       (categories: effectively_ptc + no_downstream_ejc + truncated_no_ejc).
#       Used by Panel D (stop-to-last-EJC distance density). Stop position
#       is the ref-AUG-traced comp_stop_tx_pos.
#
#   pop_ptc_plus (n=1,912 NMD / 288 Control):
#       Ref-AUG PTC+ pairs (effectively_ptc only). Used by Panels E and F.
#       Mechanism attribution combined from two sources:
#         - ref_atg_analysis$c2$attr_mechanism / attr_event for the
#           original_ptc=FALSE subset (900 NMD pairs — the ref-AUG-recovered)
#         - existing all_attr (current goal2-ptc-mechanisms output) for
#           original_ptc=TRUE pairs that are also effectively_ptc
#           (1,012 NMD pairs)

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ----
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile))
  return(normalizePath(file.path(getwd(), "data_export.R")))
})()
HERE <- dirname(script_path)
# --- DM_DIR: repointed 2026-07-26 (W10) ----------------------------------------------
# Was .P$ISOPAIR/data_mashr -- the LEGACY tree -- so this script exported panel data from
# the PUBLISHED pairs no matter what had been rebuilt. Every panel it writes is therefore
# a picture of the 1,548-pair set, which is why the committed figures are numerically
# stale and why the figure code has to move into the repo at all: the numbers moved
# (1,548 -> 1,578, n=130 -> 132), so the panels must be regenerated, not just re-pointed.
#
# Default .P$OUT/data_mashr; pass --data-dir. 01b and 02 write the same filenames there
# and only 02's pairs carry the 25% reference-share floor, so point at the floored tree.
# 2026-07-26: was .P$OUT/data_mashr, which is the UNFLOORED tree while the 25% re-scope
# is mid-flight (REFERENCE_FLOOR_PLAN R2 open). nmd_data_dir() reads .P$ISOPAIR_OUT so
# the population is a declared decision, and it prints the C4 pair count so the wrong
# tree announces itself instead of silently doubling pop_BC.
DM_DIR <- nmd_data_dir(.P)
local({
  a <- commandArgs(trailingOnly = TRUE)
  if ("--data-dir" %in% a) { i <- which(a == "--data-dir"); if (i < length(a)) DM_DIR <<- a[i + 1] }
})
DM_DIR <- normalizePath(DM_DIR, mustWork = FALSE)
if (!dir.exists(DM_DIR)) {
  stop(sprintf("DM_DIR does not exist: %s\n  Run 01 -> 01b -> 02 first, or pass --data-dir.", DM_DIR))
}
cat("data dir: ", DM_DIR, "\n", sep = "")
CACHE_DIR <- file.path(DM_DIR, "analysis_cache")
OUT_DIR <- file.path(HERE, "data")
# Rule in code: this writes, so it must not be an input root (R/load_config.R).
nmd_assert_writable(OUT_DIR, .P, "OUT_DIR")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[data_export.R] Writing to: %s\n", OUT_DIR))

# ---- Load common data ----
cat("[setup] Loading profiles + cds + ptc + ref_atg + structures\n")
profiles_c2 <- as.data.table(readRDS(file.path(DM_DIR, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM_DIR, "profiles_c4_allsamples.rds")))
cds <- as.data.table(readRDS(file.path(DM_DIR, "cds.rds")))
ptc <- as.data.table(readRDS(file.path(DM_DIR, "ptc.rds")))
structures <- as.data.table(readRDS(file.path(DM_DIR, "structures.rds")))
ref_atg <- readRDS(file.path(CACHE_DIR, "ref_atg_analysis.rds"))

# Panels D/E/F now use the all-3-ENST + coding-CDS + own-GENCODE-stop scope
# (no ref-AUG projection, no TD2 dependency). mechanism_class.R is not invoked
# by Figure 3 — Figure 4's RATIONALE.md §11 documents the per-figure scopes.
LIB <- normalizePath(file.path(HERE, "..", "..", "lib"))

# ============================================================================
# Population 1 — pop_BC: Stage 2 gene-matched (B + C)
# ============================================================================
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared_keys <- intersect(unique(make_key(profiles_c2)),
                         unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared_keys]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared_keys]
cat(sprintf("[pop_BC] Stage 2 gene-matched: NMD=%d, Control=%d\n",
            nrow(pop_BC_c2), nrow(pop_BC_c4)))

# ============================================================================
# Population 2 — all-3-ENST. Panels D/E/F universe.
#
# At all-3-ENST every isoform (reference, NMD comparator, Control comparator)
# is GENCODE-annotated, so we use each isoform's OWN GENCODE-annotated CDS
# stop position directly — no ref-AUG-projection needed. PTC status is
# determined per comparator from its own annotated stop's distance to the
# last exon-exon junction (50-nt rule).
#
# Panels B/C remain at pop_BC scope (whatever the rebuild gives); no CDS dependency so
# the all-3-ENST restriction is unnecessary for them.
# ============================================================================
is_enst <- function(x) grepl("^ENST", x)

# pop_BC + all 3 isoforms ENST + gene-matched (NMD c2 and Control c4 share gene+reference)
pop_BC_c2_E <- pop_BC_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
pop_BC_c4_E <- pop_BC_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
joint_keys <- intersect(make_key(pop_BC_c2_E), make_key(pop_BC_c4_E))
pop_BC_c2_tri <- pop_BC_c2_E[make_key(pop_BC_c2_E) %in% joint_keys]
pop_BC_c4_tri <- pop_BC_c4_E[make_key(pop_BC_c4_E) %in% joint_keys]
cat(sprintf("[all-3-ENST pre-CDS-filter] NMD=%d  Control=%d  (gene/reference combos=%d)\n",
            nrow(pop_BC_c2_tri), nrow(pop_BC_c4_tri), length(joint_keys)))

# Restrict to pairs where ALL 3 isoforms have GENCODE CDS annotation (Pete's
# text: "all three isoforms were present in GENCODE with CDS annotation").
coding_iso_ids <- cds$isoform_id[cds$coding_status == "coding"]
have_cds <- function(d) {
  d$reference_isoform_id %in% coding_iso_ids &
    d$comparator_isoform_id %in% coding_iso_ids
}
pop_BC_c2_tri_pre <- pop_BC_c2_tri  # snapshot for diagnostics
pop_BC_c2_tri <- pop_BC_c2_tri[have_cds(pop_BC_c2_tri)]
pop_BC_c4_tri <- pop_BC_c4_tri[have_cds(pop_BC_c4_tri)]
# Re-intersect on the (gene, reference) keys so both sides stay matched
joint_keys2 <- intersect(make_key(pop_BC_c2_tri), make_key(pop_BC_c4_tri))
pop_BC_c2_tri <- pop_BC_c2_tri[make_key(pop_BC_c2_tri) %in% joint_keys2]
pop_BC_c4_tri <- pop_BC_c4_tri[make_key(pop_BC_c4_tri) %in% joint_keys2]
cat(sprintf("[all-3-ENST + coding-CDS] NMD=%d  Control=%d  (gene/reference combos=%d)\n",
            nrow(pop_BC_c2_tri), nrow(pop_BC_c4_tri), length(joint_keys2)))

# ---- Compute each comparator's own GENCODE-annotated stop in tx coords ----
# cds.rds has GENCODE-annotated cds_start / cds_stop for ENST isoforms
# (SQANTI3 inherits the reference annotation for ENST-named transcripts).
suppressPackageStartupMessages({library(Isopair)})

own_stop_tx_lookup <- (function() {
  # Build a per-isoform lookup of (transcript-position-of-stop, strand)
  coding <- as.data.table(cds)[coding_status == "coding"]
  struct_idx <- match(coding$isoform_id, structures$isoform_id)
  stop_tx <- rep(NA_integer_, nrow(coding))
  for (i in seq_len(nrow(coding))) {
    si <- struct_idx[i]
    if (is.na(si)) next
    strand <- coding$strand[i]
    if (is.na(strand) || strand == "") next
    # GENCODE stop genomic position (3' end of CDS in transcript reading direction)
    stop_g <- if (strand == "+") coding$cds_stop[i] else coding$cds_start[i]
    stop_tx[i] <- Isopair::genomicToTranscript(
      stop_g, structures$exon_starts[[si]], structures$exon_ends[[si]], strand)
  }
  data.table(isoform_id = coding$isoform_id,
             own_stop_tx_pos = stop_tx,
             strand = coding$strand)
})()

# Annotate pop_BC_c2_tri / pop_BC_c4_tri with own_stop_tx_pos + last_ejc_tx_pos
augment_own_stop <- function(pop) {
  d <- merge(pop[, .(gene_id, reference_isoform_id, comparator_isoform_id)],
             own_stop_tx_lookup,
             by.x = "comparator_isoform_id", by.y = "isoform_id",
             all.x = TRUE)
  d
}
trace_c2 <- augment_own_stop(pop_BC_c2_tri)
trace_c4 <- augment_own_stop(pop_BC_c4_tri)

# ---- PTC determination using each comparator's own GENCODE stop ----
# Use structure-based junctions (.junction_positions equivalent) for the 50-nt rule.
junctions_for <- function(iso_id) {
  si <- match(iso_id, structures$isoform_id)
  if (is.na(si)) return(integer(0))
  exon_starts <- structures$exon_starts[[si]]
  exon_ends   <- structures$exon_ends[[si]]
  strand      <- structures$strand[si]
  if (is.na(strand) || strand == "") return(integer(0))
  exon_lens <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lens <- rev(exon_lens)
  if (length(exon_lens) <= 1L) return(integer(0))
  cumsum(exon_lens[-length(exon_lens)])
}

# Compute distance to last EJC and PTC flag (50 nt rule)
compute_ptc_info <- function(d) {
  d[, last_ejc_tx_pos := vapply(comparator_isoform_id,
                                 function(iso) {
                                   junc <- junctions_for(iso)
                                   if (length(junc) == 0L) NA_integer_ else max(junc)
                                 }, integer(1))]
  d[, distance := last_ejc_tx_pos - own_stop_tx_pos]
  d[, has_ptc := !is.na(distance) & distance > 50L]
  d
}
trace_c2 <- compute_ptc_info(trace_c2)
trace_c4 <- compute_ptc_info(trace_c4)

n_ptc_nmd  <- sum(trace_c2$has_ptc, na.rm = TRUE)
n_ptc_ctrl <- sum(trace_c4$has_ptc, na.rm = TRUE)
n_avail_nmd  <- sum(!is.na(trace_c2$distance))
n_avail_ctrl <- sum(!is.na(trace_c4$distance))
ptc_rate_nmd  <- 100 * n_ptc_nmd  / n_avail_nmd
ptc_rate_ctrl <- 100 * n_ptc_ctrl / n_avail_ctrl
fold_enrichment <- ptc_rate_nmd / ptc_rate_ctrl
cat(sprintf("[PTC at all-3-ENST, own GENCODE stop]\n"))
cat(sprintf("  NMD: %d/%d with PTC = %.1f%%\n", n_ptc_nmd, n_avail_nmd, ptc_rate_nmd))
cat(sprintf("  Control: %d/%d with PTC = %.1f%%\n", n_ptc_ctrl, n_avail_ctrl, ptc_rate_ctrl))
cat(sprintf("  Fold enrichment: %.1fx\n", fold_enrichment))

# pop_traceable equivalent: pairs with computable distance (Panel D)
pop_trace_c2 <- trace_c2[!is.na(distance)]
pop_trace_c4 <- trace_c4[!is.na(distance)]
# pop_ptc_plus equivalent: PTC+ pairs (Panel E/F)
pop_ptc_c2 <- trace_c2[has_ptc == TRUE]
pop_ptc_c4 <- trace_c4[has_ptc == TRUE]
# pop_ptc_c2 needs to keep the profiles_c2 metadata for attribution
# (detailed_events, n_se, etc.) — re-join from pop_BC_c2_tri
pop_ptc_c2 <- merge(pop_BC_c2_tri,
                    pop_ptc_c2[, .(comparator_isoform_id, own_stop_tx_pos,
                                   last_ejc_tx_pos, distance)],
                    by = "comparator_isoform_id")
pop_ptc_c4 <- merge(pop_BC_c4_tri,
                    pop_ptc_c4[, .(comparator_isoform_id, own_stop_tx_pos,
                                   last_ejc_tx_pos, distance)],
                    by = "comparator_isoform_id")
cat(sprintf("[pop_ptc_plus_all3ENST (GENCODE-stop-defined)]\n"))
cat(sprintf("  NMD: %d  Control: %d\n", nrow(pop_ptc_c2), nrow(pop_ptc_c4)))

# ============================================================================
# Panel B — sequence similarity NMD vs Control
#
# NOT COMPUTED HERE ANY MORE (S13, 2026-07-28). Produced by
# 05_final_report_gencode_scope_2026-07-11.Rmd :: sec1d-panelB-sequence-similarity.
# A reshape rather than a derivation -- but claim 4.2.2 is a RESULT asserted
# from this distribution, and it was stated in no report. The report now carries the
# median comparison and the Wilcoxon test behind that sentence.
# ============================================================================
panelB_path <- file.path(OUT_DIR, "panelB_sequence_similarity.tsv")
if (!file.exists(panelB_path)) {
  stop("panelB_sequence_similarity.tsv is missing. It is produced by ",
       "analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd ",
       "(chunk sec1d-panelB-sequence-similarity); render that report first. ",
       "This script no longer computes it (S13).")
}
panelB <- fread(panelB_path)
cat(sprintf("[Panel B] read %d rows from the report, not recomputed\n", nrow(panelB)))

# ============================================================================
# Panel C — splice event prevalence NMD vs Control
#
# NOT COMPUTED HERE ANY MORE (S13, 2026-07-27). The numbers now originate in
# 05_final_report_gencode_scope_2026-07-11.Rmd :: sec1d-panelC-event-prevalence,
# which writes this exact file. A result computed only in figure code is invisible
# to every report, cannot be checked by re-rendering the reports, and moves
# silently when the figure is restyled -- which is what this rule exists to stop.
#
# The migration was proven output-preserving: the report's chunk reproduces the
# TSV this block used to write, byte-identically (md5 7514c3d2a71c00f4d5ccc0d4c004af6d).
#
# ASSERT, DO NOT REGENERATE. If the file is absent the report has not been run, and
# that must fail loudly here rather than leave the panel silently stale -- a figure
# rendered from a missing input is exactly the failure mode S13 is about.
# ============================================================================
panelC_path <- file.path(OUT_DIR, "panelC_event_prevalence.tsv")
if (!file.exists(panelC_path)) {
  stop("panelC_event_prevalence.tsv is missing. It is produced by ",
       "analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd ",
       "(chunk sec1d-panelC-event-prevalence); render that report first. ",
       "This script no longer computes it (S13).")
}
panelC <- fread(panelC_path)
cat(sprintf("[Panel C] read %d event types from the report, not recomputed\n",
            nrow(panelC)))

# ============================================================================
# Panel D — stop-to-last-EJC distance, all-3-ENST scope
#
# NOT COMPUTED HERE ANY MORE (S13, 2026-07-28). Produced by
# 05_final_report_gencode_scope_2026-07-11.Rmd :: sec2a-panelD-export.
# `distance` is computed by compute_ptc_info() in the report's sec2-scope chunk;
# this was only reshaping it into long form.
# ============================================================================
panelD_path <- file.path(OUT_DIR, "panelD_stop_codon_distance.tsv")
if (!file.exists(panelD_path)) {
  stop("panelD_stop_codon_distance.tsv is missing. It is produced by ",
       "analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd ",
       "(chunk sec2a-panelD-export); render that report first. ",
       "This script no longer computes it (S13).")
}
panelD <- fread(panelD_path)
cat(sprintf("[Panel D] read %d rows from the report, not recomputed\n", nrow(panelD)))


# ============================================================================
# Panel E + F — PTC-causing event attribution + mechanism breakdown
# Population: pop_ptc_plus ∩ all-3-ENST (NMD effectively_ptc, all 3 isoforms ENST)
# Attribution source: ref-AUG-derived for ALL pairs.
# Control baseline: all detailed_events in pop_trace_c4 ∩ all-3-ENST
# ============================================================================
cat("[Panel E+F] PTC-causing attribution at pop_ptc_plus all-3-ENST scope\n")
source(file.path(HERE, "panel_e_compute.R"))

cat("\n[data_export.R] done.\n")
cat(sprintf("Population summary:\n"))
cat(sprintf("  pop_BC:                    NMD=%d Control=%d\n",
            nrow(pop_BC_c2), nrow(pop_BC_c4)))
cat(sprintf("  pop_BC (all-3-ENST):       NMD=%d Control=%d\n",
            nrow(pop_BC_c2_tri), nrow(pop_BC_c4_tri)))
cat(sprintf("  pop_traceable all-3-ENST:  NMD=%d Control=%d\n",
            nrow(pop_trace_c2), nrow(pop_trace_c4)))
cat(sprintf("  pop_ptc_plus all-3-ENST:   NMD=%d Control=%d\n",
            nrow(pop_ptc_c2), nrow(pop_ptc_c4)))
