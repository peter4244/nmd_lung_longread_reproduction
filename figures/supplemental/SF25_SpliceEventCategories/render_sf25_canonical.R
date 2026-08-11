#!/usr/bin/env Rscript
# SF25 — Twelve splice event categories detected by Isopair.
#
# Renders the 12-panel schematic using the canonical `plot_ucsc_pair` function
# from the Isopair pipeline's visualization module. Matches Yul's transcript-
# model display standard (UCSC-browser chevron introns, exon rectangles with
# CDS/UTR coloring, Reference-on-top / Comparator-on-bottom pair convention).
#
# Coordinates are the illustrative examples from the Rmd's Event Type Glossary
# section (`05_final_report_mashr.Rmd`, chunks `glossary-*`), reproduced verbatim.
#
# Output: figure_sf25_splice_event_categories.{pdf,png}

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Source the canonical visualization functions
# Repointed 2026-07-26 (W31 step 2): the helper ships in this repo; sourcing it from
# .P$ISOPAIR pulled code out of the legacy tree.
VIS <- file.path(.P$ROOT, "analysis", "isopair", "isopair_wrapper",
                 "visualization_functions.R")
source(VIS)

# Convenience helper. Panel title = human-readable event description. The
# x-axis coordinate values (0, 250, ...) are stripped, since they are not
# meaningful for these schematic exon-model panels.
glossary_plot <- function(ref_s, ref_e, comp_s, comp_e, title, subtitle) {
  # Build the pair-track panel directly from build_ucsc_track /
  # build_intron_arrows using tight y-centers (ref_y = 1.35,
  # comp_y = 1.00) so the two isoforms read as a visually paired unit.
  # The wrapper plot_ucsc_pair hardcodes ref_y = 2, comp_y = 1 which
  # leaves a large gap between rectangles; here we replicate the
  # essential geometry inline with a smaller center-to-center distance.
  ref_y  <- 1.50
  comp_y <- 1.00

  ref_track  <- build_ucsc_track(ref_s,  ref_e,  ref_y,  strand = "+")
  comp_track <- build_ucsc_track(comp_s, comp_e, comp_y, strand = "+")
  all_rects  <- rbind(cbind(ref_track,  track = "ref"),
                      cbind(comp_track, track = "comp"))

  all_coords <- c(ref_s, ref_e, comp_s, comp_e)
  fig_range  <- max(all_coords) - min(all_coords)
  chev_spacing  <- fig_range * 0.06
  chev_arm      <- chev_spacing * 0.4

  ref_intr  <- build_intron_arrows(ref_s,  ref_e,  ref_y,  strand = "+",
                                    chev_spacing = chev_spacing,
                                    chev_arm_width = chev_arm)
  comp_intr <- build_intron_arrows(comp_s, comp_e, comp_y, strand = "+",
                                    chev_spacing = chev_spacing,
                                    chev_arm_width = chev_arm)
  all_lines  <- rbind(ref_intr$lines,  comp_intr$lines)
  all_arrows <- rbind(ref_intr$arrows, comp_intr$arrows)

  x_min <- min(all_rects$xmin)
  x_max <- max(all_rects$xmax)
  x_mid <- (x_min + x_max) / 2
  x_buf <- max((x_max - x_min) * 0.05, 50)

  p <- ggplot()
  if (nrow(all_lines) > 0) {
    p <- p + geom_segment(data = all_lines,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = "gray40", linewidth = 0.4)
  }
  if (nrow(all_arrows) > 0) {
    p <- p + geom_segment(data = all_arrows,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = "gray40", linewidth = 0.3)
  }
  p <- p +
    geom_rect(data = all_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#BDC3C7", color = "black", linewidth = 0.3) +
    # Font sizes near the 10-pt docx-scale target — SF25 is rendered at
    # 10.4" native and scaled to 6.5" in the docx (factor 0.625). Single-
    # line subtitles at that size are wider than a 2.6"-native column, so
    # each subtitle wraps to 2 lines (caller supplies "\n"). Subtitle
    # docx pt = 5.0 × 2.845 × 0.625 ≈ 8.9 pt (at floor); Ref/Comp docx pt
    # = 4.6 × 2.845 × 0.625 ≈ 8.2 pt (schematic labels).
    annotate("text", x = x_mid, y = 2.30, label = subtitle,
             hjust = 0.5, vjust = 1, size = 5.0, fontface = "bold",
             lineheight = 0.9) +
    annotate("text", x = x_min, y = ref_y  + 0.18, label = "Reference",
             hjust = 0, vjust = 0, size = 4.6) +
    annotate("text", x = x_min, y = comp_y + 0.18, label = "Comparator",
             hjust = 0, vjust = 0, size = 4.6) +
    coord_cartesian(xlim = c(x_min - x_buf, x_max + x_buf),
                    ylim = c(0.55, 2.35), clip = "off") +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid    = element_blank(),
          legend.position = "none",
          plot.margin   = margin(1, 2, 1, 2),
          plot.title    = element_blank(),
          plot.subtitle = element_blank(),
          axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          axis.title    = element_blank())
  p
}

# 12 event definitions — coordinates verbatim from the Rmd Event Type Glossary.
p_alt_tss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(250, 400, 700), c(350, 550, 850),
  "Alt_TSS", "Alternative transcription\nstart site (Alt TSS)")

p_alt_tes <- glossary_plot(
  c(100, 400, 700), c(200, 550, 900),
  c(100, 400, 700), c(200, 550, 800),
  "Alt_TES", "Alternative transcription\nend site (Alt TES)")

p_a5ss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 700), c(200, 600, 850),
  "A5SS", "Alternative 5'\nsplice site (A5SS)")

p_a3ss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 350, 700), c(200, 550, 850),
  "A3SS", "Alternative 3'\nsplice site (A3SS)")

p_se <- glossary_plot(
  c(100, 400, 600, 800), c(200, 500, 700, 900),
  c(100, 400, 800),      c(200, 500, 900),
  "SE", "Skipped\nexon (SE)")

p_mi <- glossary_plot(
  c(100, 350, 550, 750, 950), c(200, 450, 650, 850, 1050),
  c(100, 950),                c(200, 1050),
  "Missing_Internal", "Missing internal\nexon(s) (Missing Internal)")

p_pir5 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 700), c(200, 620, 850),
  "Partial_IR_5", "Partial intron retention,\n5' side (Partial IR 5')")

p_pir3 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 640), c(200, 550, 850),
  "Partial_IR_3", "Partial intron retention,\n3' side (Partial IR 3')")

p_ir <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400),      c(200, 850),
  "IR", "Full intron\nretention (Full IR)")

p_ir_d5 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 430),      c(200, 850),
  "IR_diff_5", "IR with 5' boundary\ndifference (IR diff 5')")

p_ir_d3 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400),      c(200, 820),
  "IR_diff_3", "IR with 3' boundary\ndifference (IR diff 3')")

p_ir_d53 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 430),      c(200, 820),
  "IR_diff_5_3", "IR with 5' and 3' boundary\ndifferences (IR diff 5'+3')")

# 3 rows × 4 cols composite via patchwork.
composite <- (p_alt_tss | p_alt_tes | p_a5ss | p_a3ss) /
             (p_se      | p_mi      | p_pir5 | p_pir3) /
             (p_ir      | p_ir_d5   | p_ir_d3| p_ir_d53) +
             patchwork::plot_annotation(
               theme = theme(plot.margin = margin(6, 6, 6, 6))
             )

# Repointed 2026-07-26 (W31 step 2): was the LEGACY figures tree (figures/SupplementalFigures/, renamed `supplemental/` here).
HERE <- file.path(.P$ROOT, "figures", "supplemental", "SF25_SpliceEventCategories")
ggsave(file.path(HERE, "figure_sf25_splice_event_categories.pdf"),
       composite, width = 10.4, height = 8.4, device = cairo_pdf)
ggsave(file.path(HERE, "figure_sf25_splice_event_categories.png"),
       composite, width = 10.4, height = 8.4, dpi = 200)
message("wrote figure_sf25_splice_event_categories.{pdf,png}")
