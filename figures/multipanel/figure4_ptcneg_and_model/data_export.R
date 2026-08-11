#!/usr/bin/env Rscript
# ==============================================================================
# Figure 4 — §4 data export (4-panel structure)
#
# Layout (2×2):
#   Row 1 — Section A: all-3-ENST + coding-CDS (n=72/118/190)
#     A: 5'UTR length, 3 groups (PTC+ / PTC- / Control)
#     B: Longest 5'UTR ORF, 3 groups
#   Row 2 — Section C: ENST-ref + ref-AUG-traceable, 1:1 gene-matched re-intersected (n=1,050/113/1,166)
#     C: 5'UTR length, 3 groups (ref-AUG-anchored)
#     D: Longest 5'UTR ORF, 3 groups (ref-AUG-anchored)
#
# TD2-bias justification (TD2-vs-ref-AUG ORF length, paired Kozak, TD2-ATG
# position) is in the consolidated SupplementalFigures/SF34-SF35_TD2BiasEvidence/
# supplement on the unified n=900 occult-PTC scope.
#
# All NMD+/PTC- groups are MERGED (ORF match + ORF diff) — single PTC-
# bucket per Pete's guidance: "where we feel confident about the correctness
# of the CDS, just one NMD+/PTC- group."
#
# Pre-registered in ./RATIONALE.md.
# Run from the figure folder: Rscript data_export.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
  library(Isopair)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("--file=", args)]
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  normalizePath(getwd())
})()
OUT_DIR <- file.path(HERE, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this

# exported panels from the 2026-03-11 structure universe while the reports use the

# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.

DM <- nmd_data_dir(.P)
CACHE <- file.path(DM, "analysis_cache")
FASTA_PATH <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")

# Mechanism class helper (5-group → 4-group → 3-group merged here)
LIB <- normalizePath(file.path(HERE, "..", "..", "lib"))
source(file.path(LIB, "mechanism_class.R"))

# ─────────────────────────────────────────────────────────────────────
# Load
# ─────────────────────────────────────────────────────────────────────
ref_atg    <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
cds        <- as.data.table(readRDS(file.path(DM, "cds.rds")))
structures <- as.data.table(readRDS(file.path(DM, "structures.rds")))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
utr5_refaug <- readRDS(file.path(CACHE, "utr5_features_refaug.rds"))$isoform_features
utr5_all    <- readRDS(file.path(CACHE, "utr5_features_all.rds"))$isoform_features

is_enst <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

# ─────────────────────────────────────────────────────────────────────
# Universe construction
# ─────────────────────────────────────────────────────────────────────
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# Section A scope: all 3 ENST + coding-CDS + gene-matched
pop_BC_c2_E <- pop_BC_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
pop_BC_c4_E <- pop_BC_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
joint1 <- intersect(make_key(pop_BC_c2_E), make_key(pop_BC_c4_E))
pop_BC_c2_E <- pop_BC_c2_E[make_key(pop_BC_c2_E) %in% joint1]
pop_BC_c4_E <- pop_BC_c4_E[make_key(pop_BC_c4_E) %in% joint1]
coding_ids <- cds[coding_status == "coding", isoform_id]
have_cds <- function(d) d$reference_isoform_id %in% coding_ids &
                        d$comparator_isoform_id %in% coding_ids
pop_BC_c2_E2 <- pop_BC_c2_E[have_cds(pop_BC_c2_E)]
pop_BC_c4_E2 <- pop_BC_c4_E[have_cds(pop_BC_c4_E)]
joint2 <- intersect(make_key(pop_BC_c2_E2), make_key(pop_BC_c4_E2))
sec_A_c2 <- pop_BC_c2_E2[make_key(pop_BC_c2_E2) %in% joint2]
sec_A_c4 <- pop_BC_c4_E2[make_key(pop_BC_c4_E2) %in% joint2]
cat(sprintf("[Section A] all-3-ENST + coding: NMD=%d  Control=%d\n",
            nrow(sec_A_c2), nrow(sec_A_c4)))

# Section C scope: ENST-reference + ref-AUG-traceable, re-intersected on
# (gene_id, reference_isoform_id) AFTER the category filter so the NMD and
# Control panels share the same 1:1 gene-matched universe.
sec_C_c2_base <- pop_BC_c2[is_enst(reference_isoform_id)]
sec_C_c4_base <- pop_BC_c4[is_enst(reference_isoform_id)]
c2_full <- as.data.table(ref_atg$c2)[comparator_isoform_id %in%
                                       sec_C_c2_base$comparator_isoform_id]
c4_full <- as.data.table(ref_atg$c4)[comparator_isoform_id %in%
                                       sec_C_c4_base$comparator_isoform_id]
TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
# Apply TR filter to each side first, then re-intersect on (gene, reference).
c2_tr <- c2_full[category %in% TR]
c4_tr <- c4_full[category %in% TR]
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
c2_full <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c4_pt   <- merge(c4_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c2_full[, mech5 := mechanism_class(category, comp_orf_length, ref_orf_length)]
c2_full[, mech4 := mechanism_class_4(mech5)]
c2_full[, group3 := fcase(
  mech4 == "NMD+/PTC+", "NMD+/PTC+",
  mech4 %in% c("NMD+/PTC- ORF match", "NMD+/PTC- ORF diff"), "NMD+/PTC-",
  default = "Other"
)]
cat(sprintf("[Section C] ENST-reference + ref-AUG-traceable, re-intersected: NMD c2 = %d  Control c4 = %d\n",
            nrow(c2_full[group3 %in% c("NMD+/PTC+", "NMD+/PTC-")]),
            nrow(c4_pt)))

# ─────────────────────────────────────────────────────────────────────
# Section A classification — own GENCODE CDS, 3-group merged
# ─────────────────────────────────────────────────────────────────────
orf_len_lookup <- setNames(cds$orf_length, cds$isoform_id)
own_stop_genomic_lookup <- setNames(
  ifelse(cds$strand == "+", cds$cds_stop, cds$cds_start), cds$isoform_id)
own_atg_genomic_lookup <- setNames(
  ifelse(cds$strand == "+", cds$cds_start, cds$cds_stop), cds$isoform_id)
strand_lookup <- setNames(cds$strand, cds$isoform_id)

own_stop_tx <- rep(NA_integer_, nrow(cds))
last_ejc_tx <- setNames(
  vapply(seq_len(nrow(structures)), function(i) {
    n <- structures$n_exons[i]
    if (is.na(n) || n < 2L) return(NA_integer_)
    lens <- structures$exon_ends[[i]] - structures$exon_starts[[i]] + 1L
    if (structures$strand[i] == "-") lens <- rev(lens)
    as.integer(sum(lens[1:(n - 1L)]))
  }, integer(1)),
  structures$isoform_id)

for (i in seq_len(nrow(cds))) {
  if (cds$coding_status[i] != "coding") next
  iso <- cds$isoform_id[i]
  si <- match(iso, structures$isoform_id)
  if (is.na(si)) next
  strand <- cds$strand[i]
  if (is.na(strand) || strand == "") next
  stop_g <- if (strand == "+") cds$cds_stop[i] else cds$cds_start[i]
  own_stop_tx[i] <- Isopair::genomicToTranscript(
    stop_g, structures$exon_starts[[si]], structures$exon_ends[[si]], strand)
}
own_stop_tx_lookup <- setNames(own_stop_tx, cds$isoform_id)

classify_sec_A <- function(pop) {
  d <- copy(pop)
  d[, comp_stop_tx := own_stop_tx_lookup[comparator_isoform_id]]
  d[, last_ejc := last_ejc_tx[comparator_isoform_id]]
  d[, distance := last_ejc - comp_stop_tx]
  d[, comp_orf := orf_len_lookup[comparator_isoform_id]]
  d[, ref_orf  := orf_len_lookup[reference_isoform_id]]
  d[, has_ptc := !is.na(distance) & distance > 50L]
  d[, group3 := fcase(
    is.na(has_ptc), "Unclassifiable",
    has_ptc, "NMD+/PTC+",
    !has_ptc, "NMD+/PTC-"
  )]
  d
}
sec_A_c2_class <- classify_sec_A(sec_A_c2)
sec_A_c4_class <- classify_sec_A(sec_A_c4)
sec_A_c4_class[, group3 := "Control"]  # for Section A boxplot grouping
cat("Section A NMD classification:\n")
print(sec_A_c2_class[, .N, by = group3])

# ─────────────────────────────────────────────────────────────────────
# Statistical helpers
# ─────────────────────────────────────────────────────────────────────
hodges_lehmann <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2L || length(y) < 2L) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  wt <- tryCatch(wilcox.test(x, y, conf.int = TRUE, exact = FALSE),
                  error = function(e) NULL)
  if (is.null(wt)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  c(unname(wt$estimate), unname(wt$conf.int[1]), unname(wt$conf.int[2]),
    unname(wt$p.value))
}
cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0L || length(y) == 0L) return(NA_real_)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
}
pairwise_stats <- function(groups_list, lookup) {
  combos <- combn(names(groups_list), 2, simplify = FALSE)
  do.call(rbind, lapply(combos, function(p) {
    x <- lookup[groups_list[[p[1]]]]; y <- lookup[groups_list[[p[2]]]]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    hl <- hodges_lehmann(x, y)
    d  <- cliffs_delta(x, y)
    data.frame(group_x = p[1], group_y = p[2],
               n_x = length(x), n_y = length(y),
               median_x = median(x), median_y = median(y),
               hl_shift = hl[1], hl_ci_lo = hl[2], hl_ci_hi = hl[3],
               cliffs_delta = d, wilcox_p = hl[4],
               stringsAsFactors = FALSE)
  }))
}
group_descriptives <- function(groups_list, lookup, order_names) {
  out <- lapply(order_names, function(g) {
    v <- lookup[groups_list[[g]]]
    v <- v[!is.na(v)]
    data.frame(group = g, n = length(v),
               median = median(v),
               q25 = quantile(v, 0.25), q75 = quantile(v, 0.75),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ─────────────────────────────────────────────────────────────────────
# Panels A + B — Section A, all-3-ENST + coding (5'UTR length, longest uORF)
# Use utr5_features_all (TD2-anchored, which == GENCODE for ENST isoforms).
# ─────────────────────────────────────────────────────────────────────
u5_len_all <- setNames(as.numeric(utr5_all$utr5_length), utr5_all$isoform_id)
u5_orf_all <- setNames(as.numeric(utr5_all$longest_orf_nt), utr5_all$isoform_id)

sec_A_groups <- list(
  "NMD+/PTC+" = sec_A_c2_class[group3 == "NMD+/PTC+", comparator_isoform_id],
  "NMD+/PTC-" = sec_A_c2_class[group3 == "NMD+/PTC-", comparator_isoform_id],
  "Control"   = sec_A_c4_class$comparator_isoform_id
)

export_sec_A <- function(metric_name, lookup, slug) {
  desc <- group_descriptives(sec_A_groups, lookup,
                              c("NMD+/PTC+", "NMD+/PTC-", "Control"))
  stats <- pairwise_stats(sec_A_groups, lookup)
  # Long form for plotting
  long <- do.call(rbind, lapply(names(sec_A_groups), function(g) {
    v <- lookup[sec_A_groups[[g]]]
    data.frame(comparator_isoform_id = names(v), group = g, value = unname(v),
               stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$value) & long$value >= 0, ]
  setNames(long, c("comparator_isoform_id", "group", metric_name)) -> long
  fwrite(long, file.path(OUT_DIR, sprintf("panel%s_long.tsv", slug)), sep = "\t")
  fwrite(desc, file.path(OUT_DIR, sprintf("panel%s_descriptives.tsv", slug)),
         sep = "\t")
  fwrite(stats, file.path(OUT_DIR, sprintf("panel%s_pairwise.tsv", slug)),
         sep = "\t")
  cat(sprintf("\n[Section A] %s (all-3-ENST + coding)\n", metric_name))
  print(desc); print(stats)
}
export_sec_A("utr5_length_nt",  u5_len_all, "A_5utr_length")
export_sec_A("longest_5utr_orf_nt", u5_orf_all, "B_longest_5utr_orf")

# ─────────────────────────────────────────────────────────────────────
# TD2-bias evidence panels have been consolidated into the supplemental
# figure at figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/ on the unified
# n=900 occult-PTC scope. See that folder's data_export.R for the
# TD2-vs-ref-AUG ORF length, paired Kozak PWM, and TD2-ATG-position exports.
# ─────────────────────────────────────────────────────────────────────
# Panels E + F — Section C, ENST-reference + ref-AUG-traceable
# 3-group merged. Use utr5_features_refaug.rds (ref-AUG-anchored scan).
# ─────────────────────────────────────────────────────────────────────
u5_len_ra <- setNames(as.numeric(utr5_refaug$utr5_length), utr5_refaug$isoform_id)
u5_orf_ra <- setNames(as.numeric(utr5_refaug$longest_orf_nt), utr5_refaug$isoform_id)

sec_C_groups <- list(
  "NMD+/PTC+" = c2_full[group3 == "NMD+/PTC+", comparator_isoform_id],
  "NMD+/PTC-" = c2_full[group3 == "NMD+/PTC-", comparator_isoform_id],
  "Control"   = c4_pt$comparator_isoform_id
)

export_sec_C <- function(metric_name, lookup, slug) {
  desc <- group_descriptives(sec_C_groups, lookup,
                              c("NMD+/PTC+", "NMD+/PTC-", "Control"))
  stats <- pairwise_stats(sec_C_groups, lookup)
  long <- do.call(rbind, lapply(names(sec_C_groups), function(g) {
    v <- lookup[sec_C_groups[[g]]]
    data.frame(comparator_isoform_id = names(v), group = g, value = unname(v),
               stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$value) & long$value >= 0, ]
  setNames(long, c("comparator_isoform_id", "group", metric_name)) -> long
  fwrite(long, file.path(OUT_DIR, sprintf("panel%s_long.tsv", slug)), sep = "\t")
  fwrite(desc, file.path(OUT_DIR, sprintf("panel%s_descriptives.tsv", slug)),
         sep = "\t")
  fwrite(stats, file.path(OUT_DIR, sprintf("panel%s_pairwise.tsv", slug)),
         sep = "\t")
  cat(sprintf("\n[Section C] %s (ENST-ref + ref-AUG-anchored)\n", metric_name))
  print(desc); print(stats)
}
export_sec_C("utr5_length_nt",       u5_len_ra, "C_5utr_length")
export_sec_C("longest_5utr_orf_nt",  u5_orf_ra, "D_longest_5utr_orf")

cat("\n[data_export.R] done.\n")
