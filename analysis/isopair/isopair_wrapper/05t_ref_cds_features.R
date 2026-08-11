#!/usr/bin/env Rscript
# ==============================================================================
# 05t_ref_cds_features.R
#
# Compute reference-CDS structural features for ALL classified coding isoforms.
#
# For each gene, the "reference CDS" is the CDS of the dominant non-NMD isoform
# (highest mean DMSO CPM). For every other isoform in the gene, we:
#   1. Check if the reference ATG is exonic in the target
#   2. Trace the ORF from that ATG through the target's transcript
#   3. Count downstream EJCs from the traced stop codon
#   4. Compute 5'UTR features using the reference CDS boundaries
#   5. Compute 3'UTR length using the reference CDS stop position
#
# This provides a CDS-anchored feature set independent of the target isoform's
# own TransDecoder2 CDS prediction. The key question: does using the gene's
# reference reading frame improve PTC identification and NMD prediction?
#
# Inputs:
#   - data_mashr/expression_data.rds (CPM matrix)
#   - data_mashr/sample_metadata.rds (sample info for DMSO identification)
#   - data_mashr/nmd_classification.rds (NMD/non-NMD classification)
#   - data_mashr/cds.rds (TD2 CDS annotations)
#   - data_mashr/structures.rds (exon coordinates)
#   - data_mashr/gene_map.rds (isoform → gene mapping)
#   - SQANTI corrected FASTA
#
# Output:
#   - data_mashr/analysis_cache/ref_cds_features_all.rds
#     $features: data frame (one row per isoform_id) with ref-CDS features
#     $ref_lookup: gene_id → ref isoform_id mapping
#     $metadata: coverage stats, edge case counts
#
# Usage:
#   Rscript 05t_ref_cds_features.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
})

# REPOINTED 2026-07-27, same fix and same reason as 05s_orfik_scan.R.
#
# This did setwd(.P$ISOPAIR) and read "data_mashr/..." RELATIVE to it -- the LEGACY tree.
# It was missed alongside 05s: W44 recorded 05s as "the W11 class, missed when the others were
# done", and that undercounted. ref_cds_features_all.rds feeds export_rds.R -> data_prep.py ->
# the model, so reading it from the legacy tree is the stale-universe mechanism itself, on the
# structural-feature side rather than the ORF side.
#
# No fallback to the legacy tree: it is always populated, so a fallback could only ever
# silently serve the vintage the repoint exists to stop using.
DATA_MASHR <- .P$ISOPAIR_OUT
if (!dir.exists(DATA_MASHR))
  stop("intermediates root not found: ", DATA_MASHR,
       "\n  This script reads the DEPOSIT-REBUILT tree and does NOT fall back to the legacy one.")
cat("  data_mashr root:", DATA_MASHR, "\n")

FASTA_PATH <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")
OUTPUT_PATH <- file.path(.P$OUT, "data_mashr/analysis_cache/ref_cds_features_all.rds")
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

cat("=== Reference-CDS Features for All Coding Isoforms ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================
cat("Loading data...\n")

cpm_mat <- readRDS(file.path(DATA_MASHR, "expression_data.rds"))
sample_meta <- readRDS(file.path(DATA_MASHR, "sample_metadata.rds"))
nmd_class <- readRDS(file.path(DATA_MASHR, "nmd_classification.rds"))
cds <- readRDS(file.path(DATA_MASHR, "cds.rds"))
structures <- readRDS(file.path(DATA_MASHR, "structures.rds"))
gene_map <- readRDS(file.path(DATA_MASHR, "gene_map.rds"))

# Identify target population: all classified coding isoforms
nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

target_ids <- intersect(union(nmd_ids, non_nmd_ids), coding_ids)
cat("  Target coding isoforms:", length(target_ids), "\n")
cat("    NMD:", sum(target_ids %in% nmd_ids), "\n")
cat("    Non-NMD:", sum(target_ids %in% non_nmd_ids), "\n")

# ==============================================================================
# 2. Identify dominant non-NMD isoform per gene (DMSO CPM)
# ==============================================================================
cat("\nIdentifying dominant non-NMD isoform per gene...\n")

dmso_samples <- sample_meta$sample_id[sample_meta$treatment == "DMSO"]
non_nmd_coding <- intersect(non_nmd_ids, coding_ids)
non_nmd_in_cpm <- intersect(non_nmd_coding, rownames(cpm_mat))

cat("  Non-NMD coding isoforms in CPM matrix:", length(non_nmd_in_cpm), "\n")

# Mean DMSO CPM per non-NMD coding isoform
dmso_cpm <- cpm_mat[non_nmd_in_cpm, dmso_samples, drop = FALSE]
mean_dmso <- rowMeans(dmso_cpm)

# Map to genes
non_nmd_gene_map <- gene_map[gene_map$isoform_id %in% non_nmd_in_cpm, ]

# For each gene, pick the non-NMD coding isoform with highest mean DMSO CPM
ref_df <- data.frame(
  isoform_id = names(mean_dmso),
  gene_id = non_nmd_gene_map$gene_id[match(names(mean_dmso), non_nmd_gene_map$isoform_id)],
  mean_dmso_cpm = mean_dmso,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  slice_max(mean_dmso_cpm, n = 1, with_ties = FALSE) %>%
  ungroup()

ref_lookup <- setNames(ref_df$isoform_id, ref_df$gene_id)
cat("  Genes with a dominant non-NMD coding isoform:", length(ref_lookup), "\n")

# Map target isoforms to genes
target_gene_map <- gene_map[gene_map$isoform_id %in% target_ids, ]
target_genes <- target_gene_map$gene_id[match(target_ids, target_gene_map$isoform_id)]

# How many target isoforms have a reference isoform available?
has_ref <- target_genes %in% names(ref_lookup)
cat("  Target isoforms with reference available:", sum(has_ref),
    sprintf("(%.1f%%)\n", 100 * sum(has_ref) / length(has_ref)))
cat("  Target isoforms without reference:", sum(!has_ref), "\n")

# ==============================================================================
# 3. Build exon lookup and helper functions
# ==============================================================================
cat("\nBuilding exon lookups...\n")

structs_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i) {
    list(starts = structures$exon_starts[[i]],
         ends = structures$exon_ends[[i]],
         strand = structures$strand[i])
  }), structures$isoform_id)

# Check if a genomic position (and the full 3-nt codon) is exonic
is_codon_exonic <- function(pos, iso_id, strand) {
  ex <- structs_lookup[[iso_id]]
  if (is.null(ex) || is.na(pos)) return(FALSE)
  if (strand == "+") {
    positions <- c(pos, pos + 1L, pos + 2L)
  } else {
    positions <- c(pos, pos - 1L, pos - 2L)
  }
  all(sapply(positions, function(p) any(p >= ex$starts & p <= ex$ends)))
}

# Map genomic position to transcript-space (1-based)
genomic_to_transcript <- function(genomic_pos, exon_starts, exon_ends, strand) {
  if (strand == "-") {
    exon_starts <- rev(exon_starts)
    exon_ends <- rev(exon_ends)
  }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    if (strand == "+") {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (genomic_pos - es) + 1L)
    } else {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (ee - genomic_pos) + 1L)
    }
    cum <- cum + (ee - es + 1L)
  }
  NA_integer_
}

# Map transcript-space position back to genomic
transcript_to_genomic <- function(tx_pos, exon_starts, exon_ends, strand) {
  if (strand == "-") {
    exon_starts <- rev(exon_starts)
    exon_ends <- rev(exon_ends)
  }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    exon_len <- ee - es + 1L
    if (tx_pos <= cum + exon_len) {
      offset <- tx_pos - cum - 1L
      if (strand == "+") return(es + offset)
      else return(ee - offset)
    }
    cum <- cum + exon_len
  }
  NA_integer_
}

# Compute junction positions in transcript space
get_junctions <- function(exon_starts, exon_ends, strand) {
  n <- length(exon_starts)
  if (n <= 1) return(integer(0))
  exon_lengths <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  cumsum(exon_lengths)[-n]
}

stop_codons <- c("TAA", "TAG", "TGA")

# ==============================================================================
# 4. Extract transcript sequences
# ==============================================================================
cat("Extracting sequences...\n")

id_file <- tempfile(fileext = ".txt")
writeLines(target_ids, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH)
t0 <- proc.time()
raw <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t0)[3], 1), "seconds\n")
unlink(id_file)

split_lines <- strsplit(raw, "\t", fixed = TRUE)
seq_vec <- setNames(
  toupper(vapply(split_lines, `[`, character(1), 2)),
  vapply(split_lines, `[`, character(1), 1))
n_with_seq <- sum(target_ids %in% names(seq_vec))
cat("  Sequences loaded:", length(seq_vec),
    sprintf(" (targets with sequence: %d / %d)\n", n_with_seq, length(target_ids)))

# ASSERT WHAT THE EXTRACTION DID, NOT THAT THE CALL RETURNED (W122). system(intern = TRUE)
# raises a WARNING and never an error on non-zero exit, so a failed awk returns character(0),
# strsplit/vapply reduce that to an EMPTY seq_vec without complaint, and the reference-CDS
# feature build below then runs over no sequences at exit 0 -- writing an empty feature table
# that data_prep.py consumes as the model's input. A PARTIAL shortfall is legitimate and is
# reported rather than asserted; zero is not.
.awk_status <- attr(raw, "status")
if (!is.null(.awk_status) && .awk_status != 0)
  stop(sprintf("awk sequence extraction exited %s; no sequences were extracted from %s",
               .awk_status, FASTA_PATH))
if (length(target_ids) > 0 && n_with_seq == 0)
  stop(sprintf("awk returned no sequence for ANY of the %d target isoforms in %s",
               length(target_ids), FASTA_PATH))

# ==============================================================================
# 5. Core analysis: trace ref ATG through ALL target isoforms
# ==============================================================================
cat("\n=== Tracing reference ATG through all target isoforms ===\n")
cat("  (This may take 10-30 minutes for ~60K isoforms)\n\n")

# Pre-compute ref CDS 5' positions
ref_cds <- cds[match(ref_lookup, cds$isoform_id), ]
ref_cds$cds_5prime <- ifelse(ref_cds$strand == "+",
                              ref_cds$cds_start, ref_cds$cds_stop)
ref_cds_5prime <- setNames(ref_cds$cds_5prime, names(ref_lookup))
ref_strand <- setNames(ref_cds$strand, names(ref_lookup))

results <- data.frame(
  isoform_id = target_ids,
  gene_id = target_genes,
  ref_isoform_id = NA_character_,
  is_self_reference = FALSE,
  ref_atg_available = NA,
  ref_downstream_ejc = NA_integer_,
  ref_orf_length = NA_integer_,
  comp_orf_length_from_ref_atg = NA_integer_,
  ref_stop_tx_pos = NA_integer_,
  ref_stop_genomic = NA_integer_,
  ref_atg_genomic = NA_integer_,
  strand = NA_character_,
  category = NA_character_,
  stringsAsFactors = FALSE
)

# Pre-compute ref ORF length per reference isoform (avoids redundant walks)
cat("  Pre-computing reference ORF lengths...\n")
ref_orf_length_cache <- setNames(rep(NA_integer_, length(ref_lookup)), names(ref_lookup))
for (g_id_tmp in names(ref_lookup)) {
  ref_id_tmp <- ref_lookup[g_id_tmp]
  if (!ref_id_tmp %in% names(seq_vec)) next
  ref_atg_tmp <- ref_cds_5prime[g_id_tmp]
  strand_tmp <- ref_strand[g_id_tmp]
  if (is.na(ref_atg_tmp)) next
  ref_s_tmp <- structs_lookup[[ref_id_tmp]]
  if (is.null(ref_s_tmp)) next
  ref_tx_tmp <- genomic_to_transcript(ref_atg_tmp, ref_s_tmp$starts, ref_s_tmp$ends, strand_tmp)
  if (is.na(ref_tx_tmp)) next
  ref_seq_tmp <- seq_vec[ref_id_tmp]
  pos <- ref_tx_tmp + 3L
  while (pos + 2L <= nchar(ref_seq_tmp)) {
    if (substr(ref_seq_tmp, pos, pos + 2L) %in% stop_codons) {
      ref_orf_length_cache[g_id_tmp] <- pos - ref_tx_tmp
      break
    }
    pos <- pos + 3L
  }
}
cat("  Reference ORF lengths computed:", sum(!is.na(ref_orf_length_cache)),
    "/", length(ref_orf_length_cache), "\n")

t0 <- proc.time()

for (i in seq_len(nrow(results))) {
  iso_id <- results$isoform_id[i]
  g_id <- results$gene_id[i]

  # Check if gene has a reference

  if (is.na(g_id) || !g_id %in% names(ref_lookup)) {
    results$category[i] <- "no_ref_isoform"
    next
  }

  ref_id <- ref_lookup[g_id]
  results$ref_isoform_id[i] <- ref_id
  ref_atg <- ref_cds_5prime[g_id]
  strand <- ref_strand[g_id]
  results$strand[i] <- strand
  results$ref_atg_genomic[i] <- ref_atg

  if (is.na(ref_atg)) { results$category[i] <- "no_ref_cds"; next }

  # Is this isoform the reference itself?
  if (iso_id == ref_id) {
    results$is_self_reference[i] <- TRUE
    results$ref_atg_available[i] <- TRUE
    # For the reference isoform, ref-CDS features = TD2 features (by definition)
    # Still trace to get consistent feature values
  }

  # Check if ref ATG codon is exonic in target
  atg_exonic <- is_codon_exonic(ref_atg, iso_id, strand)
  results$ref_atg_available[i] <- atg_exonic

  if (!atg_exonic) { results$category[i] <- "ref_atg_lost"; next }
  if (!iso_id %in% names(seq_vec)) { results$category[i] <- "no_sequence"; next }

  target_s <- structs_lookup[[iso_id]]
  if (is.null(target_s)) { results$category[i] <- "no_structure"; next }

  # Map ref ATG to transcript position in target
  tx_pos <- genomic_to_transcript(ref_atg, target_s$starts, target_s$ends, strand)
  if (is.na(tx_pos)) { results$category[i] <- "mapping_failed"; next }

  target_seq <- seq_vec[iso_id]

  # Verify ATG at mapped position
  if (substr(target_seq, tx_pos, tx_pos + 2L) != "ATG") {
    results$category[i] <- "not_atg_in_target"
    next
  }

  # Look up pre-computed ref ORF length
  results$ref_orf_length[i] <- ref_orf_length_cache[g_id]

  # Walk codons from ref ATG in target to find stop
  comp_stop <- NA_integer_
  pos <- tx_pos + 3L
  while (pos + 2L <= nchar(target_seq)) {
    if (substr(target_seq, pos, pos + 2L) %in% stop_codons) {
      comp_stop <- pos
      break
    }
    pos <- pos + 3L
  }

  if (is.na(comp_stop)) { results$category[i] <- "no_stop_in_target"; next }

  comp_orf_len <- comp_stop - tx_pos
  results$comp_orf_length_from_ref_atg[i] <- comp_orf_len
  results$ref_stop_tx_pos[i] <- comp_stop

  # Count downstream EJCs
  junctions <- get_junctions(target_s$starts, target_s$ends, strand)
  stop_end <- comp_stop + 2L
  n_dejc <- sum(junctions > stop_end)
  results$ref_downstream_ejc[i] <- pmin(n_dejc, 5L)

  # Map stop codon back to genomic (last base of stop codon)
  stop_3prime_genomic <- transcript_to_genomic(
    comp_stop + 2L, target_s$starts, target_s$ends, strand)
  # Also map first base of stop codon
  stop_5prime_genomic <- transcript_to_genomic(
    comp_stop, target_s$starts, target_s$ends, strand)
  results$ref_stop_genomic[i] <- stop_3prime_genomic

  # Assign category
  ref_orf_len <- results$ref_orf_length[i]
  if (n_dejc > 0) {
    results$category[i] <- "effectively_ptc"
  } else if (!is.na(ref_orf_len) && comp_orf_len < ref_orf_len) {
    results$category[i] <- "truncated_no_ejc"
  } else {
    results$category[i] <- "no_downstream_ejc"
  }

  if (i %% 5000 == 0) {
    cat(sprintf("  %d / %d (%.0f sec)\n", i, nrow(results),
                (proc.time() - t0)[3]))
  }
}

elapsed_trace <- round((proc.time() - t0)[3], 1)
cat(sprintf("\n  Tracing completed in %.0f seconds\n", elapsed_trace))

# Category summary
cat("\n  Category distribution:\n")
tab <- sort(table(results$category), decreasing = TRUE)
for (nm in names(tab)) {
  cat(sprintf("    %s: %d (%.1f%%)\n", nm, tab[nm],
              100 * tab[nm] / nrow(results)))
}

# ==============================================================================
# 6. Compute ref-CDS 5'UTR features via scan5UtrFeatures
# ==============================================================================
cat("\n=== Computing ref-CDS 5'UTR features ===\n")

# For isoforms with ref ATG available, construct synthetic CDS metadata
# Convention: cds_start < cds_stop (genomic order)
#   + strand: cds_start = ATG, cds_stop = last base of stop codon
#   - strand: cds_start = last base of stop codon (smallest genomic), cds_stop = ATG

ref_available <- results[results$ref_atg_available == TRUE &
                           !is.na(results$ref_stop_genomic), ]
cat("  Isoforms with ref-CDS boundaries:", nrow(ref_available), "\n")

# Build synthetic CDS: map stop codon to genomic coordinates
# The stop codon's 3' end in transcript space (comp_stop + 2) was already
# mapped to genomic in the main loop (ref_stop_genomic).
# For the ATG, we have ref_atg_genomic.
synthetic_cds <- data.frame(
  isoform_id = ref_available$isoform_id,
  strand = ref_available$strand,
  coding_status = "coding",
  stringsAsFactors = FALSE
)

# Assign cds_start/cds_stop in genomic coordinate convention
synthetic_cds$cds_start <- ifelse(
  ref_available$strand == "+",
  ref_available$ref_atg_genomic,
  ref_available$ref_stop_genomic  # - strand: stop 3' end is smallest coord
)
synthetic_cds$cds_stop <- ifelse(
  ref_available$strand == "+",
  ref_available$ref_stop_genomic,
  ref_available$ref_atg_genomic  # - strand: ATG is largest coord
)

# Filter structures to these isoforms
structs_for_utr5 <- structures[structures$isoform_id %in% synthetic_cds$isoform_id, ]
# Filter sequences
seq_for_utr5 <- seq_vec[names(seq_vec) %in% synthetic_cds$isoform_id]

cat("  Structures available:", nrow(structs_for_utr5), "\n")
cat("  Sequences available:", length(seq_for_utr5), "\n")

# Call scan5UtrFeatures with synthetic CDS
cat("  Calling scan5UtrFeatures...\n")
t_scan <- proc.time()
ref_utr5_features <- scan5UtrFeatures(
  structures = structs_for_utr5,
  cds_metadata = synthetic_cds,
  sequences = seq_for_utr5,
  min_orf_nt = 9L,
  verbose = TRUE
)
elapsed_scan <- round((proc.time() - t_scan)[3], 1)
cat("  scan5UtrFeatures completed in", elapsed_scan, "seconds\n")
cat("  Features computed for", nrow(ref_utr5_features), "isoforms\n")

# Rename columns with ref_ prefix for clarity
ref_utr5_renamed <- ref_utr5_features %>%
  transmute(
    isoform_id,
    ref_utr5_length = utr5_length,
    ref_atg_count = n_atg,
    ref_atg_density = atg_density,
    ref_atg_strong_kozak = n_strong_kozak_atg,
    ref_atg_orphan_count = n_atg_without_orf,
    ref_uorf_count_overlapping = n_orfs_overlapping,
    ref_uorf_count_inframe = n_orfs_inframe,
    ref_uorf_count_outframe = n_orfs_outframe,
    ref_uorf_longest_nt = longest_orf_nt,
    ref_utr5_orf_coverage = pct_utr5_in_orfs,
    ref_stop_density = stop_density,
    ref_utr5_excluded = (utr5_length == 0L |
      (!is.na(atg_validated) & !atg_validated))
  )

cat("  Excluded (zero 5'UTR or ATG fail):",
    sum(ref_utr5_renamed$ref_utr5_excluded), "\n")

# ==============================================================================
# 7. Compute ref-CDS 3'UTR length
# ==============================================================================
cat("\n=== Computing ref-CDS 3'UTR length ===\n")

# Same approach as 05l_unified_model.R but using synthetic CDS stop position
ref_utr3 <- structs_for_utr5 %>%
  inner_join(synthetic_cds %>% select(isoform_id, cds_start, cds_stop, strand),
             by = c("isoform_id", "strand")) %>%
  rowwise() %>%
  mutate(
    ref_utr3_length = {
      starts <- exon_starts
      ends <- exon_ends
      if (strand == "+") {
        sum(pmax(0L, ends - pmax(starts, cds_stop)))
      } else {
        sum(pmax(0L, pmin(ends, cds_start) - starts))
      }
    }
  ) %>%
  ungroup() %>%
  transmute(isoform_id,
            ref_utr3_length,
            ref_log_utr3_length = log1p(ref_utr3_length))

cat("  3'UTR computed for", nrow(ref_utr3), "isoforms\n")
cat("  Median ref-CDS 3'UTR length:", median(ref_utr3$ref_utr3_length), "bp\n")

# ==============================================================================
# 8. Assemble final feature table
# ==============================================================================
cat("\n=== Assembling final feature table ===\n")

# Start with the core tracing results
features <- results %>%
  select(isoform_id, gene_id, ref_isoform_id, is_self_reference,
         ref_atg_available, ref_downstream_ejc, ref_orf_length,
         comp_orf_length_from_ref_atg, strand, category)

# Join 5'UTR features
features <- features %>%
  left_join(ref_utr5_renamed, by = "isoform_id")

# Join 3'UTR features
features <- features %>%
  left_join(ref_utr3, by = "isoform_id")

cat("  Total rows:", nrow(features), "\n")
cat("  With ref-CDS features (non-NA downstream_ejc):",
    sum(!is.na(features$ref_downstream_ejc)), "\n")
cat("  With 5'UTR features (non-excluded):",
    sum(!is.na(features$ref_utr5_excluded) & !features$ref_utr5_excluded), "\n")
cat("  With 3'UTR features:",
    sum(!is.na(features$ref_log_utr3_length)), "\n")

# Validation: self-reference isoforms should have consistent features
self_ref <- features[features$is_self_reference == TRUE, ]
cat("\n  Self-reference isoforms:", nrow(self_ref), "\n")
cat("  Self-ref: all ref_atg_available:", all(self_ref$ref_atg_available, na.rm = TRUE), "\n")

# ==============================================================================
# 9. Summary
# ==============================================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("  Population:", nrow(features), "classified coding isoforms\n")
cat("  Genes with reference isoform:", length(ref_lookup), "\n")

cat("\n  Coverage by NMD status:\n")
for (grp in c("NMD", "non-NMD")) {
  if (grp == "NMD") ids <- nmd_ids else ids <- non_nmd_ids
  sub <- features[features$isoform_id %in% ids, ]
  n_with <- sum(!is.na(sub$ref_downstream_ejc))
  cat(sprintf("    %s: %d / %d with ref features (%.1f%%)\n",
              grp, n_with, nrow(sub), 100 * n_with / nrow(sub)))
}

cat("\n  Category distribution (all isoforms):\n")
tab <- sort(table(features$category), decreasing = TRUE)
for (nm in names(tab)) {
  cat(sprintf("    %s: %d (%.1f%%)\n", nm, tab[nm],
              100 * tab[nm] / nrow(features)))
}

# ==============================================================================
# 10. Save
# ==============================================================================
output <- list(
  features = features,
  ref_lookup = ref_lookup,
  synthetic_cds = synthetic_cds,
  metadata = list(
    run_timestamp = Sys.time(),
    fasta_path = FASTA_PATH,
    n_target_isoforms = length(target_ids),
    n_genes_with_ref = length(ref_lookup),
    n_ref_atg_available = sum(features$ref_atg_available == TRUE, na.rm = TRUE),
    n_ref_atg_lost = sum(features$ref_atg_available == FALSE, na.rm = TRUE),
    n_no_ref = sum(features$category == "no_ref_isoform", na.rm = TRUE),
    n_with_utr5_features = sum(!is.na(features$ref_utr5_excluded) &
                                  !features$ref_utr5_excluded),
    n_with_utr3_features = sum(!is.na(features$ref_log_utr3_length)),
    elapsed_trace_sec = elapsed_trace,
    elapsed_scan_sec = elapsed_scan
  )
)

saveRDS(output, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
