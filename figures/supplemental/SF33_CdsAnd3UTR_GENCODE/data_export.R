#!/usr/bin/env Rscript
# Compute CDS length, 5'UTR length, 3'UTR length for NMD+/PTC+, NMD+/PTC-,
# and Control isoforms at the all-3-ENST + coding-CDS scope (130/130, 4-CT re-scope).
#
# Reproduces the filter cascade in
# /Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model/data_export.R

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
})

# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this

# exported panels from the 2026-03-11 structure universe while the reports use the

# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.

ISOPAIR_DIR <- dirname(nmd_data_dir(.P))
DM <- file.path(ISOPAIR_DIR, "data_mashr")

dt_c2  <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
dt_c4  <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
cds    <- readRDS(file.path(DM, "cds.rds"))
struct <- readRDS(file.path(DM, "structures.rds"))
ref_atg <- readRDS(file.path(DM, "analysis_cache/ref_atg_analysis.rds"))
ra_c2 <- as.data.table(ref_atg$c2)[, .(comparator_isoform_id,
                                       utr3_to_tx_end_nt,
                                       utr3_via_non_ptc_stop_nt)]
ra_c4 <- as.data.table(ref_atg$c4)[, .(comparator_isoform_id,
                                       utr3_to_tx_end_nt,
                                       utr3_via_non_ptc_stop_nt)]

# pop_BC: pairs whose (gene_id, reference_isoform_id) appears in both c2 and c4
pop_BC_keys <- merge(
  unique(dt_c2[, .(gene_id, reference_isoform_id)]),
  unique(dt_c4[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
pop_BC_c2 <- merge(dt_c2, pop_BC_keys, by = c("gene_id", "reference_isoform_id"))
pop_BC_c4 <- merge(dt_c4, pop_BC_keys, by = c("gene_id", "reference_isoform_id"))

# all-3-ENST per side
ent <- function(x) startsWith(x, "ENST")
sub_c2 <- pop_BC_c2[ent(reference_isoform_id) & ent(comparator_isoform_id)]
sub_c4 <- pop_BC_c4[ent(reference_isoform_id) & ent(comparator_isoform_id)]
# joint (all-3-ENST gene-matched)
key_all3 <- merge(
  unique(sub_c2[, .(gene_id, reference_isoform_id)]),
  unique(sub_c4[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
j1_c2 <- merge(sub_c2, key_all3, by = c("gene_id", "reference_isoform_id"))
j1_c4 <- merge(sub_c4, key_all3, by = c("gene_id", "reference_isoform_id"))

cat(sprintf("After all-3-ENST: c2=%d, c4=%d, unique gene/ref combos=%d\n",
            nrow(j1_c2), nrow(j1_c4), nrow(key_all3)))

# Coding-everywhere filter
cds_dt <- as.data.table(cds)
coding_ids <- cds_dt[coding_status == "coding", isoform_id]
keep_c2 <- j1_c2[(reference_isoform_id %in% coding_ids) &
                 (comparator_isoform_id %in% coding_ids)]
keep_c4 <- j1_c4[(reference_isoform_id %in% coding_ids) &
                 (comparator_isoform_id %in% coding_ids)]
# Re-intersect on gene+reference
key_keep <- merge(
  unique(keep_c2[, .(gene_id, reference_isoform_id)]),
  unique(keep_c4[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
final_c2 <- merge(keep_c2, key_keep, by = c("gene_id", "reference_isoform_id"))
final_c4 <- merge(keep_c4, key_keep, by = c("gene_id", "reference_isoform_id"))

cat(sprintf("After coding + re-intersect: c2=%d, c4=%d (expect 130 each)\n",
            nrow(final_c2), nrow(final_c4)))

# Compute own-GENCODE-stop PTC classification on the COMPARATOR
# Need: own_stop_tx_pos and last_ejc_tx_pos for each comparator
struct_dt <- as.data.table(struct)
genomic_to_tx <- function(gp, es, ee, s) {
  if (s == "-") { es <- rev(es); ee <- rev(ee) }
  cum <- 0L
  for (j in seq_along(es)) {
    if (s == "+" && gp >= es[j] && gp <= ee[j]) return(cum + (gp - es[j]) + 1L)
    if (s == "-" && gp >= es[j] && gp <= ee[j]) return(cum + (ee[j] - gp) + 1L)
    cum <- cum + (ee[j] - es[j] + 1L)
  }
  NA_integer_
}

compute_features <- function(comp_ids) {
  out <- data.table(comparator_isoform_id = comp_ids)
  c_cds <- cds_dt[match(comp_ids, isoform_id)]
  c_str <- struct_dt[match(comp_ids, isoform_id)]

  # Own GENCODE-stop genomic position
  own_stop_genomic <- ifelse(c_cds$strand == "+", c_cds$cds_stop, c_cds$cds_start)
  own_atg_genomic  <- ifelse(c_cds$strand == "+", c_cds$cds_start, c_cds$cds_stop)

  own_stop_tx <- integer(nrow(out))
  own_atg_tx  <- integer(nrow(out))
  last_ejc_tx <- integer(nrow(out))
  tx_len      <- integer(nrow(out))

  for (i in seq_len(nrow(out))) {
    if (is.na(own_stop_genomic[i])) {
      own_stop_tx[i] <- NA_integer_; own_atg_tx[i] <- NA_integer_
      last_ejc_tx[i] <- NA_integer_; tx_len[i] <- NA_integer_
      next
    }
    es <- c_str$exon_starts[[i]]; ee <- c_str$exon_ends[[i]]
    s  <- c_str$strand[i];        n_ex <- length(es)
    exon_lens <- ee - es + 1L
    tx_len[i] <- sum(exon_lens)
    # last EJC = end of penultimate exon in transcript order
    if (n_ex >= 2L) {
      if (s == "+") last_ejc_tx[i] <- sum(exon_lens[1:(n_ex - 1L)])
      else          last_ejc_tx[i] <- sum(rev(exon_lens)[1:(n_ex - 1L)])
    } else {
      last_ejc_tx[i] <- NA_integer_   # single-exon = no EJC
    }
    own_stop_tx[i] <- genomic_to_tx(own_stop_genomic[i], es, ee, s)
    own_atg_tx[i]  <- genomic_to_tx(own_atg_genomic[i],  es, ee, s)
  }

  out[, `:=`(tx_len = tx_len,
             own_atg_tx = own_atg_tx,
             own_stop_tx = own_stop_tx,
             last_ejc_tx = last_ejc_tx)]
  # 5'UTR = ATG position - 1
  out[, utr5_nt := pmax(0L, own_atg_tx - 1L)]
  # CDS length = stop - start + 1 (+3 for stop codon, conventionally included)
  out[, cds_nt := pmax(0L, own_stop_tx - own_atg_tx + 1L)]
  # 3'UTR = tx_len - stop_tx_pos
  out[, utr3_nt := pmax(0L, tx_len - own_stop_tx)]
  # PTC+: own stop > 50 nt past last EJC (last_ejc_tx - own_stop_tx > 50)
  # Convention: positive distance = stop UPSTREAM of last EJC (PTC direction)
  out[, distance_to_last_ejc := last_ejc_tx - own_stop_tx]
  out[, is_ptc_plus := !is.na(distance_to_last_ejc) & distance_to_last_ejc > 50L]
  out
}

nmd_feat  <- compute_features(final_c2$comparator_isoform_id)
ctrl_feat <- compute_features(final_c4$comparator_isoform_id)

# Translation-based 3'UTR = tx_len - own_GENCODE_stop_tx_pos (consistent with
# own-GENCODE-stop PTC classification used in this scope).
# Already computed above as `utr3_nt`.
nmd_feat[,  utr3_to_tx_end_nt := utr3_nt]
ctrl_feat[, utr3_to_tx_end_nt := utr3_nt]

# Non-PTC-stop 3'UTR:
#   * For PTC- and Control (own GENCODE stop is non-PTC by classification): equal
#     to the translation-based measure (the own stop IS non-PTC).
#   * For PTC+ (own GENCODE stop is a PTC): walk past the PTC to find the next
#     non-PTC stop. Use ref_atg's pre-computed non_ptc_stop_tx_pos where available.
#     (At this all-3-ENST + coding-CDS scope, the comparator's own GENCODE ATG
#     and the ref-AUG generally coincide for same-gene isoforms, so the
#     ref_atg-based walk is a close proxy.)
ra_c2_meta <- as.data.table(ref_atg$c2)[, .(comparator_isoform_id,
                                            non_ptc_stop_tx_pos,
                                            utr3_via_non_ptc_stop_nt)]
nmd_feat <- merge(nmd_feat, ra_c2_meta, by = "comparator_isoform_id", all.x = TRUE)
# For PTC- (own stop non-PTC): non-PTC-stop measure = translation-based
nmd_feat[is_ptc_plus == FALSE, utr3_via_non_ptc_stop_nt := utr3_nt]
# For PTC+ where ref_atg field is missing, leave NA
ctrl_feat[, utr3_via_non_ptc_stop_nt := utr3_nt]

nmd_feat[,  group := ifelse(is_ptc_plus, "NMD+/PTC+", "NMD+/PTC-")]
ctrl_feat[, group := "Control"]

all_feat <- rbind(nmd_feat, ctrl_feat, fill = TRUE)
all_feat[, group := factor(group, levels = c("NMD+/PTC+", "NMD+/PTC-", "Control"))]

cat(sprintf("\nGroup counts: PTC+=%d (NMD), PTC-=%d (NMD), Control=%d\n",
            nmd_feat[is_ptc_plus == TRUE, .N],
            nmd_feat[is_ptc_plus == FALSE, .N],
            ctrl_feat[, .N]))

# Per-group descriptives for the three quantities
describe <- function(x) sprintf("median=%.0f (IQR %.0f-%.0f)",
                                median(x, na.rm = TRUE),
                                quantile(x, 0.25, na.rm = TRUE),
                                quantile(x, 0.75, na.rm = TRUE))
cat("\n=== CDS length (nt) ===\n")
cat(sprintf("  NMD+/PTC+: %s\n", describe(nmd_feat[is_ptc_plus == TRUE,  cds_nt])))
cat(sprintf("  NMD+/PTC-: %s\n", describe(nmd_feat[is_ptc_plus == FALSE, cds_nt])))
cat(sprintf("  Control:   %s\n", describe(ctrl_feat$cds_nt)))

cat("\n=== 5'UTR length (nt) ===\n")
cat(sprintf("  NMD+/PTC+: %s\n", describe(nmd_feat[is_ptc_plus == TRUE,  utr5_nt])))
cat(sprintf("  NMD+/PTC-: %s\n", describe(nmd_feat[is_ptc_plus == FALSE, utr5_nt])))
cat(sprintf("  Control:   %s\n", describe(ctrl_feat$utr5_nt)))

cat("\n=== 3'UTR length, translation-based (own stop -> tx end) ===\n")
cat(sprintf("  NMD+/PTC+: %s\n", describe(nmd_feat[is_ptc_plus == TRUE,  utr3_to_tx_end_nt])))
cat(sprintf("  NMD+/PTC-: %s\n", describe(nmd_feat[is_ptc_plus == FALSE, utr3_to_tx_end_nt])))
cat(sprintf("  Control:   %s\n", describe(ctrl_feat$utr3_to_tx_end_nt)))

cat("\n=== 3'UTR length, non-PTC-stop based (walk past PTCs) ===\n")
cat(sprintf("  NMD+/PTC+: %s  (non-NA: %d/%d)\n",
            describe(nmd_feat[is_ptc_plus == TRUE,  utr3_via_non_ptc_stop_nt]),
            sum(!is.na(nmd_feat[is_ptc_plus == TRUE,  utr3_via_non_ptc_stop_nt])),
            nmd_feat[is_ptc_plus == TRUE, .N]))
cat(sprintf("  NMD+/PTC-: %s  (non-NA: %d/%d)\n",
            describe(nmd_feat[is_ptc_plus == FALSE, utr3_via_non_ptc_stop_nt]),
            sum(!is.na(nmd_feat[is_ptc_plus == FALSE, utr3_via_non_ptc_stop_nt])),
            nmd_feat[is_ptc_plus == FALSE, .N]))
cat(sprintf("  Control:   %s  (non-NA: %d/%d)\n",
            describe(ctrl_feat$utr3_via_non_ptc_stop_nt),
            sum(!is.na(ctrl_feat$utr3_via_non_ptc_stop_nt)),
            nrow(ctrl_feat)))

# Pairwise Wilcoxon tests
wtest <- function(a, b, lab) {
  w <- wilcox.test(a, b, conf.int = TRUE)
  hl <- as.numeric(w$estimate)
  cat(sprintf("  %s: HL shift=%.0f, p=%.2e\n", lab, hl, w$p.value))
}

cat("\n=== Wilcoxon comparisons ===\n")
cat("[CDS length]\n")
wtest(nmd_feat[is_ptc_plus == FALSE, cds_nt],
      nmd_feat[is_ptc_plus == TRUE,  cds_nt], "PTC- vs PTC+")
wtest(nmd_feat[is_ptc_plus == FALSE, cds_nt], ctrl_feat$cds_nt, "PTC- vs Control")
wtest(nmd_feat[is_ptc_plus == TRUE,  cds_nt], ctrl_feat$cds_nt, "PTC+ vs Control")

cat("[5'UTR length]\n")
wtest(nmd_feat[is_ptc_plus == FALSE, utr5_nt],
      nmd_feat[is_ptc_plus == TRUE,  utr5_nt], "PTC- vs PTC+")
wtest(nmd_feat[is_ptc_plus == FALSE, utr5_nt], ctrl_feat$utr5_nt, "PTC- vs Control")
wtest(nmd_feat[is_ptc_plus == TRUE,  utr5_nt], ctrl_feat$utr5_nt, "PTC+ vs Control")

cat("[3'UTR length, translation-based]\n")
wtest(nmd_feat[is_ptc_plus == FALSE, utr3_to_tx_end_nt],
      nmd_feat[is_ptc_plus == TRUE,  utr3_to_tx_end_nt], "PTC- vs PTC+")
wtest(nmd_feat[is_ptc_plus == FALSE, utr3_to_tx_end_nt], ctrl_feat$utr3_to_tx_end_nt, "PTC- vs Control")
wtest(nmd_feat[is_ptc_plus == TRUE,  utr3_to_tx_end_nt], ctrl_feat$utr3_to_tx_end_nt, "PTC+ vs Control")

cat("[3'UTR length, non-PTC-stop based]\n")
wtest(nmd_feat[is_ptc_plus == FALSE, utr3_via_non_ptc_stop_nt],
      nmd_feat[is_ptc_plus == TRUE,  utr3_via_non_ptc_stop_nt], "PTC- vs PTC+")
wtest(nmd_feat[is_ptc_plus == FALSE, utr3_via_non_ptc_stop_nt], ctrl_feat$utr3_via_non_ptc_stop_nt, "PTC- vs Control")
wtest(nmd_feat[is_ptc_plus == TRUE,  utr3_via_non_ptc_stop_nt], ctrl_feat$utr3_via_non_ptc_stop_nt, "PTC+ vs Control")

# Save for further use
# Repointed 2026-07-26 (W31 step 2): was the LEGACY figures tree (figures/SupplementalFigures/, renamed `supplemental/` here).
SUPP_DIR <- file.path(.P$ROOT, "figures", "supplemental", "SF33_CdsAnd3UTR_GENCODE")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SUPP_DIR, "data"), showWarnings = FALSE)
fwrite(all_feat[, .(comparator_isoform_id, group, cds_nt, utr5_nt,
                    utr3_to_tx_end_nt, utr3_via_non_ptc_stop_nt)],
       file.path(SUPP_DIR, "data", "features_190_subset_long.tsv"), sep = "\t")
cat(sprintf("\nSaved per-isoform feature table -> %s/data/features_190_subset_long.tsv\n", SUPP_DIR))

# Per-measure descriptives + pairwise tests
emit_per_measure <- function(measure_col, slug) {
  groups <- c("NMD+/PTC+", "NMD+/PTC-", "Control")
  desc <- rbindlist(lapply(groups, function(g) {
    vals <- all_feat[group == g, get(measure_col)]
    vals <- vals[!is.na(vals)]
    data.table(group = g, n = length(vals),
               median = median(vals),
               q25 = quantile(vals, 0.25),
               q75 = quantile(vals, 0.75))
  }))
  fwrite(desc, file.path(SUPP_DIR, "data",
                        sprintf("%s_descriptives.tsv", slug)), sep = "\t")
  pairs <- list(c("NMD+/PTC-", "NMD+/PTC+"),
                c("NMD+/PTC-", "Control"),
                c("NMD+/PTC+", "Control"))
  pw <- rbindlist(lapply(pairs, function(p) {
    a <- all_feat[group == p[1], get(measure_col)]; a <- a[!is.na(a)]
    b <- all_feat[group == p[2], get(measure_col)]; b <- b[!is.na(b)]
    w <- suppressWarnings(wilcox.test(a, b, conf.int = TRUE))
    cliffs_n <- length(a) * length(b)
    cliffs <- (sum(outer(a, b, ">")) - sum(outer(a, b, "<"))) / cliffs_n
    data.table(group_x = p[1], group_y = p[2],
               n_x = length(a), n_y = length(b),
               median_x = median(a), median_y = median(b),
               hl_shift = as.numeric(w$estimate),
               hl_ci_lo = as.numeric(w$conf.int[1]),
               hl_ci_hi = as.numeric(w$conf.int[2]),
               cliffs_delta = cliffs,
               wilcox_p = w$p.value)
  }))
  fwrite(pw, file.path(SUPP_DIR, "data",
                      sprintf("%s_pairwise.tsv", slug)), sep = "\t")
  invisible(NULL)
}
emit_per_measure("cds_nt", "panelA_cds_length")
emit_per_measure("utr5_nt", "panelB_utr5_length")
emit_per_measure("utr3_to_tx_end_nt", "panelC_utr3_translation")
emit_per_measure("utr3_via_non_ptc_stop_nt", "panelD_utr3_non_ptc_stop")
cat(sprintf("Saved per-panel descriptives + pairwise TSVs in %s/data/\n", SUPP_DIR))

# ============================================================================
# Section C (ref-AUG-projected CDS / 3'UTR) — READ, do not compute.
#
# S13: no result may originate in `figures/`.
#
# Until 2026-07-28 this block rebuilt the Section C cohort ITSELF — it repeated the
# sec3-scope construction (pop_BC -> ENST-reference -> ref-AUG-traceable category ->
# re-intersect) and re-ran mechanism_class() -> group3, then emitted the feature table and
# the six panel TSVs. The canonical report performs the SAME construction in its own
# `sec3-scope` chunk, so one population was defined in two places and only one of them
# could be checked by re-rendering a report (S4). Claim 4.7.5's operative half — that the
# 3'UTR effect is limited to the ref-AUG-traceable set — was produced entirely here,
# outside the report layer.
#
# They now originate in
#   analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd
#   :: sec3d-refaug-cds-utr3-export
# and this script only reads what that chunk writes.
#
# PROVEN OUTPUT-PRESERVING BEFORE THE MOVE. Both constructions were run against the same
# tree and compared directly: identical 833/833 cohort with an identical key set, an
# identical 1,666-row feature table, and identical descriptives AND pairwise tables for all
# three measures (cds_nt, utr3_to_tx_end_nt, utr3_via_non_ptc_stop_nt).
#
# THE `n=1,166` IN THE OLD HEADER WAS STALE. It is a legacy-tree vintage, not a different
# scope; the deposit-native cohort is 833 per arm. The filename
# `features_1166_subset_long.tsv` preserves that stale number and is kept only so the
# SF36 reader does not break — it is misleading and should be renamed when SF36 is next
# touched.
#
# ASSERT, DO NOT REGENERATE. A missing file means the report has not been run, and that
# must stop the pipeline rather than leave a panel silently stale.
# ============================================================================
.secC_files <- c("features_1166_subset_long.tsv",
                 "panelD_cds_length_refaug_descriptives.tsv",
                 "panelD_cds_length_refaug_pairwise.tsv",
                 "panelE_utr3_translation_refaug_descriptives.tsv",
                 "panelE_utr3_translation_refaug_pairwise.tsv",
                 "panelF_utr3_non_ptc_stop_refaug_descriptives.tsv",
                 "panelF_utr3_non_ptc_stop_refaug_pairwise.tsv")

for (.f in .secC_files) {
  .path <- file.path(SUPP_DIR, "data", .f)
  if (!file.exists(.path)) {
    stop(sprintf(paste0(
      "%s is missing. It is produced by analysis/isopair/isopair_wrapper/",
      "05_final_report_gencode_scope_2026-07-11.Rmd (chunk sec3d-refaug-cds-utr3-export); ",
      "render that report first. This script no longer computes it (S13)."), .f))
  }
}
all_feat_C <- fread(file.path(SUPP_DIR, "data", "features_1166_subset_long.tsv"))
cat(sprintf("[Section C] read %d rows from the report, not recomputed\n", nrow(all_feat_C)))
cat("  group sizes:\n"); print(all_feat_C[, .N, by = group])
