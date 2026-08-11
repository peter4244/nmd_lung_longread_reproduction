#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5 — n=1,166 Subset 2 (ref-AUG-traceable) isoform list + Path B-strict
# uORF attention computation.
#
# This is the broader pairwise scope from Figure 4 Section C — gene-matched
# 1:1 (NMD c2 vs Control c4), ENST-reference, ref-AUG-traceable (i.e.,
# mechanism_class ∈ {effectively_ptc, no_downstream_ejc, truncated_no_ejc}).
# Total: 1,166 NMD comparator pairs + 1,166 Control comparator pairs.
#
# Subgroup structure at this scope:
#   NMD+/PTC+         = category == "effectively_ptc"  (~1,050)
#   NMD+/PTC- retained = category ∈ {no_downstream_ejc, truncated_no_ejc}  (~116)
#   Control            = all 1,166 c4 comparators
#   (NMD+/PTC- lost not in this scope by construction — ref AUG is required
#   to project, and "lost" means it doesn't project)
#
# Joins to results_4ct/uorf_attention_metrics.tsv (v1) and applies
# Path B-strict structural reclassification locally using the priority-slot
# data (uorf_features_in_priority_slots.tsv) where available.
#
# Outputs:
#   data/subset2_refaug_isoforms.tsv       (one row per comparator with group)
#   data/panelG_uorf_attention_refaug.tsv  (per-isoform v1 + Path B-strict)
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0)
    return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
OUT_DIR <- file.path(HERE, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this

# exported panels from the 2026-03-11 structure universe while the reports use the

# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.

DM <- nmd_data_dir(.P)
CACHE <- file.path(DM, "analysis_cache")
LIB <- normalizePath(file.path(HERE, "..", "..", "lib"))
source(file.path(LIB, "mechanism_class.R"))

# MODEL_RESULTS, not file.path(.P$MODEL, "results_4ct"). This is the ONE necessary file with
# an unguarded read into the model repo, which is why the rerun's isolation spec exempts
# .P$MODEL from the roots that must not resolve -- and that exemption admitted the whole
# PUBLISHED results_4ct to this script. Every §5 step taking --results-dir is invoked with
# results_4ct_dn (D6), so this one silently exported published values into a Figure 5 panel
# and would have "agreed" with the paper because it read it. W76.
#
# A paths.yml KEY rather than a bespoke env var, because paths.yml is the single resolver
# (NMD_<KEY> > paths.yml > repo-relative) and a second mechanism is how roots drift apart.
# Override for a deposit-native run: NMD_MODEL_RESULTS=../NMD_orf_model_v5_4ct/results_4ct_dn
MODEL_RES <- .P$MODEL_RESULTS
cat(sprintf("[model-results] %s\n", MODEL_RES))

# tx_summary.tsv is a STAGE C feature table (D105), not a model artifact, and with stage E
# skipped it is NOT beside the model outputs: MODEL_RESULTS is the deposit's trained
# artifacts and the feature tables are in stage C's results dir. FEATURES defaults to
# MODEL_RESULTS so the normal flow is unchanged; the deposit-native harness sets NMD_FEATURES.
# NO FALLBACK BETWEEN THE TWO: a search that tries one root then the other would silently take
# whichever happened to exist and report success, which is the failure the runbook names.
FEATURES <- .P$FEATURES
cat(sprintf("[features]      %s\n", FEATURES))

# ── Load ──
cat("Loading Isopair profiles + ref_atg cache...\n")
ref_atg     <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

is_enst  <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

# ── Build pop_BC (gene-matched) ──
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# ── Section C scope: ENST-ref + ref-AUG-traceable, re-intersected 1:1 ──
sec_C_c2_base <- pop_BC_c2[is_enst(reference_isoform_id)]
sec_C_c4_base <- pop_BC_c4[is_enst(reference_isoform_id)]
c2_full <- as.data.table(ref_atg$c2)[comparator_isoform_id %in% sec_C_c2_base$comparator_isoform_id]
c4_full <- as.data.table(ref_atg$c4)[comparator_isoform_id %in% sec_C_c4_base$comparator_isoform_id]
TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
c2_tr <- c2_full[category %in% TR]
c4_tr <- c4_full[category %in% TR]
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
c2_refaug <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c4_refaug <- merge(c4_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))

# 4-CT re-scope + 25% ref-share floor (2026-07-11): n=1,166 -> n=819 in the PUBLISHED run.
#
# THE COUNT IS REPORTED, NOT ASSERTED (2026-08-03). This was
# `stopifnot(nrow(c2_refaug) == 819L, nrow(c4_refaug) == 819L)` -- the PUBLISHED number, welded
# into a script whose job is to rebuild published numbers. A clean-room run from the deposit
# produces 833/833, which is not an error: the claim ledger already records it, as 4.7.4's
# "768/833 = 92.2% vs published 756/819 = 92.3%". So the assertion could only ever pass on the
# published data and halted every correct regeneration by construction.
#
# WHAT THE CHECK WAS ACTUALLY PROTECTING is kept, because the hazard is real: a subset that
# silently loses rows, or two arms that stop matching, would corrupt every downstream group
# comparison without any error. Those still halt. Only the pin to one vintage is gone.
PUBLISHED_N <- 819L   # the published 4-CT count, for comparison only -- never a gate
cat(sprintf("[Subset 2] NMD c2 = %d, Control c4 = %d  (published run: %d/%d)\n",
            nrow(c2_refaug), nrow(c4_refaug), PUBLISHED_N, PUBLISHED_N))
if (nrow(c2_refaug) != PUBLISHED_N || nrow(c4_refaug) != PUBLISHED_N) {
  cat(sprintf("[Subset 2] NOTE: differs from the published count by %+d/%+d. Expected for a\n",
              nrow(c2_refaug) - PUBLISHED_N, nrow(c4_refaug) - PUBLISHED_N),
      "           deposit-native rebuild; report the count with the figure rather than",
      "assuming\n            either number.\n", sep = "")
}
# The arms must match each other and be non-empty. These are the failures that would corrupt the
# group comparison silently, so they remain fatal.
stopifnot("subset 2 arms disagree -- the paired comparison would be meaningless" =
            nrow(c2_refaug) == nrow(c4_refaug),
          "subset 2 is empty -- the ref-AUG join produced nothing" =
            nrow(c2_refaug) > 0L)

# ── Assign group labels ──
c2_refaug[, group := fifelse(category == "effectively_ptc",
                             "NMD+/PTC+", "NMD+/PTC- retained")]
c2_refaug[, arm := "NMD"]
c4_refaug[, group := "Control"]
c4_refaug[, arm := "Control"]

out <- rbind(
  c2_refaug[, .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, category)],
  c4_refaug[, .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, category)]
)

cat("\nSubgroup counts:\n")
print(out[, .N, by = group])

fwrite(out, file.path(OUT_DIR, "subset2_refaug_isoforms.tsv"), sep = "\t")

# ── Join with model uORF metrics + priority-slot file ──
mtx <- fread(file.path(MODEL_RES, "uorf_attention_metrics.tsv"),
             select = c("isoform_id", "uorf_attention_frac", "top_is_uorf",
                        "main_orf_rank", "split", "label", "prob"))
slot <- fread(file.path(MODEL_RES, "uorf_features_in_priority_slots.tsv"))
tx_file <- file.path(FEATURES, "tx_summary.tsv")
if (!file.exists(tx_file))
  stop("stage C feature table not found: ", tx_file,
       "\n  tx_summary.tsv is a stage C output (D105), not a deposit file.",
       "\n  set NMD_FEATURES to the directory export_rds.R wrote with --results-dir.")
tx   <- fread(tx_file, select = c("isoform_id", "tx_length"))

# Path B-strict per slot
slot <- merge(slot, tx, by = "isoform_id")
slot[, orf_start_approx := frac_position * tx_length]
slot[, frac_stop := (orf_start_approx + orf_length) / tx_length]
slot[, is_strict := orf_length < 200 &
                     frac_position < 0.5 &
                     frac_stop < 0.6 &
                     kozak_score >= 1]

# Per isoform: sum attention over Path B-strict slots
# Isoforms NOT in slot file → strict attention = 0 (no flagged-uORF attention by v1)
pb <- slot[, .(strict_attn = sum(attention[is_strict])), by = isoform_id]

# Join: n=1,166 list + v1 metrics + strict attention
panG_refaug <- merge(out, mtx, by.x = "comparator_isoform_id", by.y = "isoform_id",
                    all.x = TRUE)
panG_refaug <- merge(panG_refaug, pb, by.x = "comparator_isoform_id", by.y = "isoform_id",
                    all.x = TRUE)
panG_refaug[is.na(strict_attn), strict_attn := 0]

cat(sprintf("\nModel coverage: %d / %d (%.1f%%) have v1 metrics from model\n",
            sum(!is.na(panG_refaug$uorf_attention_frac)),
            nrow(panG_refaug),
            100 * mean(!is.na(panG_refaug$uorf_attention_frac))))

# ── Subgroup table ──
gord <- c("NMD+/PTC+", "NMD+/PTC- retained", "Control")
have <- panG_refaug[!is.na(uorf_attention_frac)]

cat("\n=== Path B-strict vs v1: uORF attention by subgroup at n=1,166 ===\n\n")
tab <- have[, .(
  n              = .N,
  v1_mean        = round(mean(uorf_attention_frac), 3),
  v1_median      = round(median(uorf_attention_frac), 3),
  strict_mean    = round(mean(strict_attn), 3),
  strict_median  = round(median(strict_attn), 3),
  v1_top_is_uorf_pct     = round(100 * mean(top_is_uorf, na.rm = TRUE), 1),
  strict_any_attn_pct    = round(100 * mean(strict_attn > 0), 1),
  strict_attn_gt05_pct   = round(100 * mean(strict_attn > 0.05), 1),
  strict_attn_gt20_pct   = round(100 * mean(strict_attn > 0.20), 1)
), by = group][order(factor(group, levels = gord))]
print(tab)

cat("\n=== Pairwise Wilcoxon of Path B-strict attention ===\n")
pairs <- combn(gord, 2, simplify = FALSE)
pw <- rbindlist(lapply(pairs, function(p) {
  x <- have[group == p[1], strict_attn]
  y <- have[group == p[2], strict_attn]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.table(group_x = p[1], group_y = p[2],
             nx = length(x), ny = length(y),
             wilcox_p = signif(wt$p.value, 3))
}))
print(pw)

# Descriptives — for the panel labels
desc <- have[, .(n = .N,
                 median = round(median(strict_attn), 4),
                 q25 = round(quantile(strict_attn, 0.25), 4),
                 q75 = round(quantile(strict_attn, 0.75), 4),
                 mean = round(mean(strict_attn), 4),
                 pct_any = round(100 * mean(strict_attn > 0), 1),
                 pct_gt05 = round(100 * mean(strict_attn > 0.05), 1),
                 pct_gt20 = round(100 * mean(strict_attn > 0.20), 1)),
             by = group][order(factor(group, levels = gord))]

fwrite(panG_refaug, file.path(OUT_DIR, "panelG_uorf_attention_refaug.tsv"), sep = "\t")
fwrite(pw,         file.path(OUT_DIR, "panelG_uorf_attention_refaug_pairwise.tsv"), sep = "\t")
fwrite(desc,       file.path(OUT_DIR, "panelG_uorf_attention_refaug_descriptives.tsv"), sep = "\t")

cat("\nWrote:\n  ", file.path(OUT_DIR, "subset2_refaug_isoforms.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_refaug.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_refaug_pairwise.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_refaug_descriptives.tsv"), "\n",
    sep="")
