#!/usr/bin/env Rscript
# =============================================================================
# 01_extract_our_isoforms.R
#
# Build the per-isoform input table for the NMD predictor comparison.
# Output schema (columns):
#   gene_id, comparator_isoform_id, reference_isoform_id, subclass,
#   in_last_exon, distance_to_start_nt, exon_length_at_stop_nt,
#   within_50nt_of_last_ej, mashr_posterior_mean_logfc, our_model_prob, chr
#
# - Cohort: n = 1,166 ref-AUG-traceable scope per METHODS.md §2.1
#     1,050 NMD+/PTC+ (c2 effectively_ptc) +
#       116 NMD+/PTC- (c2 no_downstream_ejc + truncated_no_ejc) +
#     1,166 Control comparators (c4, gene-matched on (gene, reference))
# - Features: the four NMDetective-B rules, computed from each comparator's
#     long-read-observed transcript structure (not canonical lookup).
# - Gold standard: per-isoform mashr posterior mean log fold-change under
#     SMG1i (4-CT scope), aggregated across the four manuscript cell types
#     by max (matches the union convention of the `nmd_responsive` flag).
# - Our model predictions: pulled from predictions_atg500_stop500.tsv.
# =============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
setwd(HERE)   # scripts still resolve siblings relative to themselves

DM    <- file.path(.P$OUT, "data_mashr")
CACHE <- file.path(DM, "analysis_cache")
MASHR <- file.path(file.path(.P$OUT, "mashr_isoform"), "mashr_isoform_posterior_means_2026.3.10.csv")
# Full-cohort predictions (split="all"): emitted by run_infer_all.py on Explorer
# (scoring train + val + test + test_paralog isoforms in the model's H5).
# Columns: isoform_id, chr, h5_split, label, logit, prob.
# REPOINTED 2026-08-06 onto the deposited 1000nt model. The old name is not in the deposit at
# all, so this chain -- and SF43 behind it -- could not run.
PREDS <- file.path(.P$DEPOSIT, "model", "predictions_atg1000_stop1000_seed42_all.tsv")
# Output root. 01 wrote to PC_OUT without defining it -- so this script has never run to
# completion in this form. Defined here to match 02/03/04.
PC_OUT <- file.path(.P$OUT, "predictor_comparison")
dir.create(PC_OUT, showWarnings = FALSE, recursive = TRUE)
nmd_assert_writable(PC_OUT, .P, "PC_OUT")

# 2026.8.6, not 2026.7.11: these outputs are the DEPOSITED 1000nt model rescored.
# Leaving the July stamp would put new data in a file named for the old vintage.
DATESTAMP <- "2026.8.6"

# ── Load ──
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
ref_atg     <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
cds         <- as.data.table(readRDS(file.path(DM, "cds.rds")))
structures  <- as.data.table(readRDS(file.path(DM, "structures.rds")))

# ── Build the n=1,166 ref-AUG-traceable cohort (matches figure4 data_export.R) ──
is_enst  <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# Restrict to ENST-reference + ref-AUG-traceable categories on the c2 side
sec_C_c2_base <- pop_c2[is_enst(reference_isoform_id)]
sec_C_c4_base <- pop_c4[is_enst(reference_isoform_id)]

c2_full <- as.data.table(ref_atg$c2)[comparator_isoform_id %in%
                                       sec_C_c2_base$comparator_isoform_id]
c4_full <- as.data.table(ref_atg$c4)[comparator_isoform_id %in%
                                       sec_C_c4_base$comparator_isoform_id]

TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
c2_tr <- c2_full[category %in% TR]
c4_tr <- c4_full[category %in% TR]  # match Fig 4 construction (n=1,166 scope)

# Re-intersect on (gene_id, reference_isoform_id) after category filter
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id")
)
c2_cohort <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c4_cohort <- merge(c4_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))

cat(sprintf("Cohort: NMD comparators = %d, Control comparators = %d\n",
            nrow(c2_cohort), nrow(c4_cohort)))

# ── Subclass assignment ──
c2_cohort[, subclass := fcase(
  category == "effectively_ptc",                            "NMD+/PTC+",
  category %in% c("no_downstream_ejc", "truncated_no_ejc"), "NMD+/PTC-",
  default = NA_character_
)]
c4_cohort[, subclass := "Control"]

# Print subclass counts
cat("\nSubclass counts:\n")
print(c2_cohort[, .N, by = subclass])
print(c4_cohort[, .N, by = subclass])

# ── Per-isoform feature computation: the four NMDetective-B rules ──
# Stop-codon transcript position chosen per subclass:
#   NMD+/PTC+:   ref-AUG-projected stop  (comp_stop_tx_pos)
#   NMD+/PTC-:   ref-AUG-projected stop  (same column)
#   Control:     natural stop from cds.rds (own GENCODE-annotated stop in
#                comparator-transcript coordinates, derived below)

# Pre-compute lookups from structures.rds
tx_len_lookup <- setNames(
  vapply(seq_len(nrow(structures)), function(i) {
    s <- structures$exon_starts[[i]]; e <- structures$exon_ends[[i]]
    if (is.null(s) || length(s) == 0L) return(NA_integer_)
    as.integer(sum(e - s + 1L))
  }, integer(1)),
  structures$isoform_id
)

# Per-isoform last-exon transcript-coordinate boundary and last-EJ position
build_struct_lookups <- function() {
  lst <- list()
  for (i in seq_len(nrow(structures))) {
    iso  <- structures$isoform_id[i]
    s    <- structures$exon_starts[[i]]
    e    <- structures$exon_ends[[i]]
    if (is.null(s) || length(s) == 0L) next
    lens <- e - s + 1L
    if (structures$strand[i] == "-") lens <- rev(lens)
    cs   <- cumsum(lens)       # cumulative tx position at end of each exon
    n    <- length(lens)
    lst[[iso]] <- list(
      n_exons         = n,
      last_exon_start = if (n >= 1) c(0, cs[-n])[n] + 1L else NA_integer_,
      last_ej_tx_pos  = if (n >= 2) as.integer(cs[n - 1L]) else NA_integer_,
      exon_end_tx     = as.integer(cs)
    )
  }
  lst
}
cat("\nBuilding per-isoform structure lookups ...\n")
struct_lkp <- build_struct_lookups()

# Comparator-side natural-stop transcript position for Control comparators
own_stop_tx_lookup <- setNames(
  vapply(seq_len(nrow(cds)), function(i) {
    if (cds$coding_status[i] != "coding") return(NA_integer_)
    iso <- cds$isoform_id[i]
    si  <- match(iso, structures$isoform_id)
    if (is.na(si)) return(NA_integer_)
    strand <- cds$strand[i]
    if (is.na(strand) || strand == "") return(NA_integer_)
    stop_g <- if (strand == "+") cds$cds_stop[i] else cds$cds_start[i]
    Isopair::genomicToTranscript(
      stop_g, structures$exon_starts[[si]], structures$exon_ends[[si]], strand)
  }, integer(1)),
  cds$isoform_id
)

# Function to compute the 4 features from a (isoform_id, stop_tx_pos) pair
compute_features <- function(iso, stop_tx) {
  if (is.na(stop_tx) || is.null(struct_lkp[[iso]])) {
    return(list(in_last_exon = NA, distance_to_start_nt = NA_integer_,
                exon_length_at_stop_nt = NA_integer_,
                within_50nt_of_last_ej = NA))
  }
  s <- struct_lkp[[iso]]
  # in_last_exon: stop_tx > start of last exon
  in_last <- stop_tx >= s$last_exon_start
  # distance to start: tx coord of start codon to stop codon
  #   We use stop_tx itself as "distance from tx start", since
  #   `comp_stop_tx_pos` is measured from the start of the spliced mRNA
  #   (i.e., the AUG of the projected ORF). For Controls, the natural
  #   stop's tx position is also relative to tx start. For NMDetective's
  #   "DistanceToStart" we want distance from AUG to stop; we approximate
  #   this as `stop_tx - tx-position-of-AUG`. For c2 isoforms,
  #   ref_atg_exonic_in_comp is the AUG tx position; we'll pass that in
  #   from the caller. Here `stop_tx` is treated as "distance from AUG to
  #   stop" directly per the input convention.
  dist_to_start <- as.integer(stop_tx)
  # exon containing the stop: first exon-end >= stop_tx
  ee <- s$exon_end_tx
  ei <- which(ee >= stop_tx)[1]
  if (is.na(ei) || ei == 0L) {
    exon_len <- NA_integer_
  } else {
    exon_start <- if (ei == 1L) 1L else ee[ei - 1L] + 1L
    exon_len   <- as.integer(ee[ei] - exon_start + 1L)
  }
  # within 50 nt of last EJ: last_ej - stop_tx, but with the canonical
  # convention "PTC if last_ej - stop_tx > 50". So "within 50 nt" means
  # last_ej - stop_tx <= 50 (i.e., not a PTC by this rule).
  if (is.na(s$last_ej_tx_pos)) {
    within_50nt <- TRUE  # single-exon transcripts: stop is by definition not before any EJ
  } else {
    within_50nt <- (s$last_ej_tx_pos - stop_tx) <= 50L
  }
  list(in_last_exon = in_last,
       distance_to_start_nt = dist_to_start,
       exon_length_at_stop_nt = exon_len,
       within_50nt_of_last_ej = within_50nt)
}

# Apply per-comparator
make_feature_table <- function(d, stop_col) {
  out <- vector("list", nrow(d))
  for (i in seq_len(nrow(d))) {
    iso  <- d$comparator_isoform_id[i]
    stp  <- as.integer(d[[stop_col]][i])
    f    <- compute_features(iso, stp)
    out[[i]] <- data.table(
      comparator_isoform_id  = iso,
      stop_tx_pos_used       = stp,
      in_last_exon           = f$in_last_exon,
      distance_to_start_nt   = f$distance_to_start_nt,
      exon_length_at_stop_nt = f$exon_length_at_stop_nt,
      within_50nt_of_last_ej = f$within_50nt_of_last_ej
    )
  }
  rbindlist(out)
}

cat("\nComputing features for NMD comparators ...\n")
nmd_feats <- make_feature_table(c2_cohort, "comp_stop_tx_pos")

# Controls: their "stop_tx_pos" is the natural-stop tx position
c4_cohort[, comp_natural_stop_tx := own_stop_tx_lookup[comparator_isoform_id]]
cat(sprintf("Controls with natural-stop tx position: %d of %d\n",
            sum(!is.na(c4_cohort$comp_natural_stop_tx)), nrow(c4_cohort)))
cat("Computing features for Controls ...\n")
ctrl_feats <- make_feature_table(c4_cohort, "comp_natural_stop_tx")

# ── Stitch the cohort + subclass + features together ──
nmd_rows <- merge(c2_cohort[, .(gene_id, reference_isoform_id,
                                comparator_isoform_id, subclass)],
                   nmd_feats, by = "comparator_isoform_id")
ctrl_rows <- merge(c4_cohort[, .(gene_id, reference_isoform_id,
                                 comparator_isoform_id, subclass)],
                    ctrl_feats, by = "comparator_isoform_id")
cohort <- rbind(nmd_rows, ctrl_rows)

# ── Gold standard: mashr posterior mean log fold-change ──
# 4-CT scope. Aggregation across CTs: MAX (matches the union convention
# behind `nmd_responsive == TRUE`).
mashr_pm <- fread(MASHR)
# Resolved rather than hard-coded -- see R/mashr_columns.R (silent-zero failure mode).
ct_cols <- unname(nmd_mash_cols(c("AT2","LAE","FB","MV"), mashr_pm))
mashr_pm[, mashr_posterior_mean_logfc := pmax(get(ct_cols[1]),
                                               get(ct_cols[2]),
                                               get(ct_cols[3]),
                                               get(ct_cols[4]),
                                               na.rm = TRUE)]
cohort <- merge(cohort,
                 mashr_pm[, .(comparator_isoform_id = txid,
                              mashr_posterior_mean_logfc)],
                 by = "comparator_isoform_id", all.x = TRUE)

# ── Our model probability ──
preds <- fread(PREDS, select = c("isoform_id", "chr", "h5_split", "prob"))
cohort <- merge(cohort,
                 preds[, .(comparator_isoform_id = isoform_id, chr,
                            h5_split,
                            our_model_prob = prob)],
                 by = "comparator_isoform_id", all.x = TRUE)

cat(sprintf("\nFinal cohort: %d rows\n", nrow(cohort)))
cat(sprintf("  with mashr log-FC:    %d\n", sum(!is.na(cohort$mashr_posterior_mean_logfc))))
cat(sprintf("  with our_model_prob:  %d\n", sum(!is.na(cohort$our_model_prob))))

# ── Write ──
out_file <- file.path(PC_OUT, sprintf("isoforms_%s.tsv", DATESTAMP))
fwrite(cohort, out_file, sep = "\t")
cat(sprintf("\nWrote %s (%d rows × %d cols)\n",
            basename(out_file), nrow(cohort), ncol(cohort)))

cat("\nSubclass × mashr-log-FC summary:\n")
print(cohort[, .(n = .N,
                  median_logfc = median(mashr_posterior_mean_logfc, na.rm = TRUE),
                  q25 = quantile(mashr_posterior_mean_logfc, 0.25, na.rm = TRUE),
                  q75 = quantile(mashr_posterior_mean_logfc, 0.75, na.rm = TRUE)),
              by = subclass])
