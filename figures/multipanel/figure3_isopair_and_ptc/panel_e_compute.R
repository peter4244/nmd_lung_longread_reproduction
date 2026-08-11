# Figure 3 Panels E and F — READ, do not compute.
#
# Sourced by data_export.R. S13: no result may originate in `figures/`.
#
# Until 2026-07-27 this file computed both panels itself: it read fc_c2_allsamples.rds and
# fw_c2_allsamples.rds, re-ran attribute_ptc_events() / attribute_3utr_splice(), rebuilt the
# Control event baseline, and wrote panelE_ptc_event_attribution.tsv and
# panelF_mechanism_breakdown.tsv. The canonical report performed the SAME attribution in its
# own `sec2b-attribution` chunk from the same two caches — so the numbers behind Figure 3 E/F
# existed in two places, and only one of them could be checked by re-rendering a report.
#
# They now originate in
#   analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd
#   :: sec2b-panelEF-export
# and this file only reads what that chunk writes.
#
# PROVEN OUTPUT-PRESERVING — both panels regenerate byte-identically after the migration:
#   panelE 8ccd62cf7c39673b673dd69dc1d8e06d
#   panelF f06288f9c2aad5f645aa2cbc19e3659f
#
# THE FIRST ATTEMPT DID NOT MATCH, and the byte-identity test is what caught it: the
# IR_diff -> IR merge must be applied to `ptc_causing_event` at the SOURCE, before
# event_short is derived, so that Panel E's proportions and Panel F's mechanism breakdown
# both see it. Collapsing only the Control side left Panel E with IR n_ptc=3 plus a stray
# IR_diff_3 row, against the figure's IR n_ptc=4. Reading the code had not revealed that;
# comparing the bytes did.
#
# ASSERT, DO NOT REGENERATE. A missing file means the report has not been run, and that must
# stop the pipeline rather than leave a panel silently stale — a figure drawn from a missing
# or out-of-date input is precisely the failure S13 exists to prevent.

.panelEF <- c(panelE = "panelE_ptc_event_attribution.tsv",
              panelF = "panelF_mechanism_breakdown.tsv")

for (.nm in names(.panelEF)) {
  .path <- file.path(OUT_DIR, .panelEF[[.nm]])
  if (!file.exists(.path)) {
    stop(sprintf(paste0(
      "%s is missing. It is produced by analysis/isopair/isopair_wrapper/",
      "05_final_report_gencode_scope_2026-07-11.Rmd (chunk sec2b-panelEF-export); ",
      "render that report first. This script no longer computes it (S13)."),
      .panelEF[[.nm]]))
  }
}

panelE <- fread(file.path(OUT_DIR, .panelEF[["panelE"]]))
panelF <- fread(file.path(OUT_DIR, .panelEF[["panelF"]]))

cat(sprintf("  [E] read %d event types from the report, not recomputed\n",
            nrow(panelE)))
cat(sprintf("  [F] read %d (event_type x mechanism) rows; total attributed=%d\n",
            nrow(panelF), sum(panelF$n)))
