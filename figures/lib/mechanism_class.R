# mechanism_class.R — derived 5-group classification for §4 / Figure 4.
#
# Pre-registered in figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md §3.
# Decoupled from ref_atg_analysis.rds (NOT a cache field) so the classification
# can evolve without invalidating Figure 3 caches.
#
# Usage:
#   source("figures/lib/mechanism_class.R")
#   c2$mechanism_class <- mechanism_class(c2$category,
#                                          c2$comp_orf_length, c2$ref_orf_length)
#   c2$ref_source <- ref_source(c2$reference_isoform_id)
#
# Inputs:
#   category         — chr, raw Isopair::traceReferenceAtg output:
#                      "effectively_ptc","no_downstream_ejc","truncated_no_ejc",
#                      "ref_atg_lost","no_ref_cds","mapping_failed", or other.
#   comp_orf_length  — int, length (nt) of the ref-AUG-projected ORF in the comparator.
#   ref_orf_length   — int, length (nt) of the ref AUG ORF in the reference.
#
# Output: chr factor with levels (in canonical display order):
#   "NMD+/PTC+", "NMD+/PTC- ORF same", "NMD+/PTC- ORF extended",
#   "NMD+/PTC- ORF truncated", "Ref AUG absent"
#
# Design decision: the classification is 100% ref-AUG-derived. We do NOT use
# Isopair's `original_ptc` field, because that field is TD2-derived from the
# comparator's CDS, and TD2 has documented PTC-avoidance bias (Figure 3 finding).
# A small set of pairs (≈13 at the ENST-only c2 scope) where TD2 called PTC+ but
# ref-AUG-tracing did not reproduce the PTC determination are routed by their
# ref-AUG category (ndj/trc) and ORF-length comparison, not by TD2's call.
#
# Notes:
#   - "Ref AUG absent" pools ref_atg_lost, no_ref_cds, mapping_failed, and any
#     unrecognized category. Excluded from §4 mechanism inference.
#   - This function is deterministic, vectorized, and has no external dependencies
#     beyond base R.

MECHANISM_CLASS_LEVELS <- c(
  "NMD+/PTC+",
  "NMD+/PTC- ORF same",
  "NMD+/PTC- ORF extended",
  "NMD+/PTC- ORF truncated",
  "Ref AUG absent"
)

# 4-group MANUSCRIPT scheme — used for §4 mechanism comparison.
# Collapses "ORF extended" + "ORF truncated" into "ORF diff".
# Pre-registered in figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md §3.
MECHANISM_CLASS_4_LEVELS <- c(
  "NMD+/PTC+",
  "NMD+/PTC- ORF match",
  "NMD+/PTC- ORF diff",
  "Ref AUG absent"
)

REF_AUG_PRESENT_CATEGORIES <- c(
  "effectively_ptc", "no_downstream_ejc", "truncated_no_ejc"
)

mechanism_class <- function(category, comp_orf_length, ref_orf_length) {
  stopifnot(length(category) == length(comp_orf_length),
            length(category) == length(ref_orf_length))

  # ORF-length difference (positive = comparator longer)
  orf_diff <- comp_orf_length - ref_orf_length

  # Group assignment (mutually exclusive, in evaluation order)
  cls <- rep(NA_character_, length(category))

  cls[is.na(category) | !category %in% REF_AUG_PRESENT_CATEGORIES] <-
    "Ref AUG absent"

  # PTC+ — purely ref-AUG-derived. Does NOT include TD2's `original_ptc` field.
  cls[!is.na(category) & category == "effectively_ptc"] <- "NMD+/PTC+"

  # PTC- (ref-AUG) subgroups — ALL ndj/trc pairs routed by orf_diff regardless
  # of TD2's call. TD2 PTC+ on ndj/trc → still ref-AUG-PTC- under §4 framework.
  ndj_mask <- !is.na(category) & category == "no_downstream_ejc"
  trc_mask <- !is.na(category) & category == "truncated_no_ejc"

  cls[ndj_mask & !is.na(orf_diff) & orf_diff == 0] <- "NMD+/PTC- ORF same"
  cls[ndj_mask & !is.na(orf_diff) & orf_diff > 0]  <- "NMD+/PTC- ORF extended"
  cls[trc_mask & !is.na(orf_diff) & orf_diff < 0]  <- "NMD+/PTC- ORF truncated"

  # ndj/trc pairs with NA orf_diff (typically because ref_orf_length is NA in
  # the cache — the reference has no annotated ORF length to compare against)
  # still have a known ref-AUG-traceable category (ndj/trc → PTC- by
  # definition). The ORF-comparison sub-detail is unknown, so we route
  # ndj-NA to "ORF same" and trc-NA to "ORF truncated"; both collapse to
  # "NMD+/PTC-" under mechanism_class_4(). This corrects a 2026-06-15
  # fall-through that previously sent these pairs to "Ref AUG absent",
  # which contradicted their (definitionally ref-AUG-traceable) category.
  cls[ndj_mask & is.na(orf_diff)] <- "NMD+/PTC- ORF same"
  cls[trc_mask & is.na(orf_diff)] <- "NMD+/PTC- ORF truncated"

  # Anything still NA falls to "Ref AUG absent" (unexpected category values
  # only). Defensive default.
  cls[is.na(cls)] <- "Ref AUG absent"

  factor(cls, levels = MECHANISM_CLASS_LEVELS)
}

# ─────────────────────────────────────────────────────────────────────
# mechanism_class_4: 4-group MANUSCRIPT scheme.
# Collapses "ORF extended" + "ORF truncated" → "ORF diff".
# All other groups pass through. Use this for §4 mechanism comparison;
# the 5-group classification is retained as internal Panel-D landscape detail.
# ─────────────────────────────────────────────────────────────────────
mechanism_class_4 <- function(cls) {
  # cls is the output of mechanism_class() (factor with 5 levels)
  cls_char <- as.character(cls)
  out <- ifelse(cls_char %in% c("NMD+/PTC- ORF extended",
                                "NMD+/PTC- ORF truncated"),
                "NMD+/PTC- ORF diff",
                ifelse(cls_char == "NMD+/PTC- ORF same",
                       "NMD+/PTC- ORF match",
                       cls_char))
  factor(out, levels = MECHANISM_CLASS_4_LEVELS)
}

# ─────────────────────────────────────────────────────────────────────
# ref_source: classify reference isoform annotation source.
# Used to restrict §4 primary analyses to GENCODE-annotated references.
# ─────────────────────────────────────────────────────────────────────
ref_source <- function(reference_isoform_id) {
  out <- ifelse(grepl("^ENST", reference_isoform_id), "GENCODE", "TD2")
  factor(out, levels = c("GENCODE", "TD2"))
}

# ─────────────────────────────────────────────────────────────────────
# tx_len_delta_concordant: TES-concordant flag for 3'UTR primary analysis.
# Pre-registered threshold: |delta| <= 50 nt (RATIONALE.md §4.2).
# ─────────────────────────────────────────────────────────────────────
TES_CONCORDANT_THRESHOLD_NT <- 50L

tes_concordant <- function(comp_tx_len, ref_tx_len,
                           threshold = TES_CONCORDANT_THRESHOLD_NT) {
  abs(comp_tx_len - ref_tx_len) <= threshold
}
