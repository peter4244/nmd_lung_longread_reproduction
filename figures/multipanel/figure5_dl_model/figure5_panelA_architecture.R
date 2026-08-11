#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5 Panel A — Model architecture schematic.
#
# Sources the upstream architecture builder in the linked model repo
# (~/claude_projects/NMD_orf_model_v5_4ct/make_architecture_figure.R), which
# builds a ggplot `p` containing the v5_4ct architecture diagram (branch
# colors: teal=ATG, coral=Stop, gold=Structural; warm grey for fusion /
# attention / classification). Saves PDF + PNG into this folder.
#
# The upstream script is a static hand-drawn schematic (no data inputs),
# so we treat it as a one-shot render rather than porting to matplotlib.
#
# Run from this folder:  Rscript figure5_panelA_architecture.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(ggplot2)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0)
    return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()

# VENDORED 2026-08-10. This read file.path(.P$MODEL, ...) and so drew Panel A only on a machine
# with the model repo checked out beside this one -- clean-room job 9049007 died here with
# "file.exists(UPSTREAM) is not TRUE". A citable repository cannot depend on a sibling checkout.
# The copy carries the upstream commit and sha256 in its header, and tools/check_vendored.py fails
# when it drifts from the recorded hash, so a divergence is loud rather than silent.
UPSTREAM <- file.path(HERE, "make_architecture_figure.R")
stopifnot(file.exists(UPSTREAM))

# Source — `p` is left in the global environment after running.
source(UPSTREAM, local = FALSE)

out_pdf <- file.path(HERE, "figure5_panelA_architecture.pdf")
out_png <- file.path(HERE, "figure5_panelA_architecture.png")

# Source script uses coord_fixed(ratio=1) over xlim c(-3,18), ylim c(5.5,30)
# → coordinate aspect 21:24.5 ≈ 0.857. ggplot font sizes are in mm and don't
# auto-scale, so at small output canvases the per-element text overruns the
# boxes. Render at 14 x 16.33 in so the text has room (aspect ≈ 0.857).
ggsave(out_pdf, plot = p, width = 14, height = 16.33, units = "in", device = cairo_pdf)
ggsave(out_png, plot = p, width = 14, height = 16.33, units = "in", dpi = 300)

cat(sprintf("Saved:\n  %s\n  %s\n", out_pdf, out_png))
