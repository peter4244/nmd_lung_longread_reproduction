#!/usr/bin/env Rscript
# SF30 — GAIN/LOSS direction by splice-event type, NMD vs Control (post 25% floor).
#
# Standalone export replacing the DEPRECATED legacy-Rmd producer. Reproduces the
# legacy `goal1-table1c-gain-loss` logic on the FLOORED, pop_BC-matched profiles:
# each pair's $detailed_events carries event_type + direction (GAIN/LOSS from the
# comparator's perspective). Count by event_type x direction for NMD (c2) and
# Control (c4), % GAIN, Fisher on the 2x2 (GAIN/LOSS x NMD/Control) per event.
# See REFERENCE_FLOOR_SF_INVENTORY.md.
#
# Output (SF30_GainDirectionByEvent/data/): table1c_gain_loss_events.csv.

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressMessages({library(data.table); library(dplyr); library(tidyr)})
# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this
# exported panels from the 2026-03-11 structure universe while the reports use the
# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.
DM <- nmd_data_dir(.P)
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)))),
          error = function(e) ".")
OUT <- file.path(HERE, "data"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

p2 <- readRDS(file.path(DM, "profiles_c2_allsamples.rds"))
p4 <- readRDS(file.path(DM, "profiles_c4_allsamples.rds"))
mk <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared <- intersect(unique(mk(p2)), unique(mk(p4)))
p2 <- p2[mk(p2) %in% shared, ]; p4 <- p4[mk(p4) %in% shared, ]   # pop_BC-matched

get_gain_loss <- function(profiles, label) {
  ev <- do.call(rbind, profiles$detailed_events)
  ev |> group_by(event_type, direction) |> summarise(n = n(), .groups = "drop") |>
    group_by(event_type) |> mutate(total = sum(n), pct = round(100 * n / total, 1)) |>
    ungroup() |> mutate(comparison = label)
}
gl <- rbind(get_gain_loss(p2, "NMD"), get_gain_loss(p4, "Control"))
wide <- gl |> pivot_wider(id_cols = event_type,
                          names_from = c(comparison, direction),
                          values_from = c(n, pct), values_fill = 0) |> arrange(event_type)

cols <- c("event_type", "n_NMD_GAIN", "n_NMD_LOSS", "pct_NMD_GAIN",
          "n_Control_GAIN", "n_Control_LOSS", "pct_Control_GAIN")
disp <- as.data.frame(wide)[, intersect(cols, names(wide)), drop = FALSE]
disp$fisher_p <- NA_real_; disp$direction_diff <- NA_character_
for (i in seq_len(nrow(disp))) {
  ng <- disp$n_NMD_GAIN[i]; nl <- disp$n_NMD_LOSS[i]
  cg <- disp$n_Control_GAIN[i]; cl <- disp$n_Control_LOSS[i]
  if (ng + nl > 0 && cg + cl > 0) {
    ft <- tryCatch(fisher.test(matrix(c(ng, nl, cg, cl), nrow = 2)), error = function(e) NULL)
    if (!is.null(ft)) {
      disp$fisher_p[i] <- signif(ft$p.value, 3)
      disp$direction_diff[i] <- if (ft$p.value > 0.05) "Non-significant"
        else if (100*ng/(ng+nl) > 100*cg/(cg+cl)) "More GAIN in NMD" else "More LOSS in NMD"
    }
  }
}
write.csv(disp, file.path(OUT, "table1c_gain_loss_events.csv"), row.names = FALSE)
cat(sprintf("SF30: %d event types over %d NMD / %d Control pop_BC pairs\n",
            nrow(disp), nrow(p2), nrow(p4)))
