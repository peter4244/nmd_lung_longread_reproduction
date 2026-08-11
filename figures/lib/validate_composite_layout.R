#!/usr/bin/env Rscript
# validate_composite_layout.R
#
# Validates spatial layout of patchwork composite figures.
#
# Checks:
#   1. Tag-title collision: panel tag overlaps a title at the top of the
#      same panel (especially inside wrap_elements)
#   2. Tag-tag overlap: two tags share the same space
#   3. Missing or duplicate panel tags
#   4. Tag font consistency across panels
#
# Strategy: the gtable layout gives us reliable cell positions. Tags live
# in `tag-N` cells. Titles either live in their own `title-N` cells
# (safe — patchwork allocates space) or inside panel body / `full-N` cells
# (potential collision if positioned at the top). We check for collisions
# by examining cell adjacency and text npc positions within cells.
#
# Usage:
#   source("figures/lib/validate_composite_layout.R")
#   combined <- (pA | pB) / (pC | pD) +
#     plot_annotation(tag_levels = "A")
#   result <- validate_composite_layout(combined, fig_width = 16, fig_height = 18)
#
# Requires: ggplot2, patchwork, grid

# ---------------------------------------------------------------------------
# Helper: recursively find text grobs in a grob/gtable tree.
# Returns a list of (label, fontsize, fontface, y_npc, path).
# ---------------------------------------------------------------------------
.find_text_grobs <- function(grob, path = "") {
  results <- list()

  if (inherits(grob, "text")) {
    for (i in seq_along(grob$label)) {
      lab <- grob$label[i]
      if (is.null(lab) || is.na(lab) || !nzchar(trimws(lab))) next
      results[[length(results) + 1]] <- list(
        label    = lab,
        fontsize = if (!is.null(grob$gp$fontsize)) grob$gp$fontsize[1] else NA_real_,
        fontface = if (!is.null(grob$gp$fontface)) grob$gp$fontface[1] else NA,
        y_npc    = as.numeric(grob$y[i]),
        path     = path
      )
    }
  }

  children <- NULL
  if (inherits(grob, "gtable") && !is.null(grob$grobs)) {
    children <- grob$grobs
  } else if (inherits(grob, "gTree") && !is.null(grob$children)) {
    children <- grob$children
  }
  if (!is.null(children)) {
    for (i in seq_along(children)) {
      nm <- if (!is.null(names(children))) names(children)[i] else as.character(i)
      results <- c(results, .find_text_grobs(children[[i]], paste0(path, "/", nm)))
    }
  }

  results
}

# ---------------------------------------------------------------------------
# Helper: classify a text grob found inside a panel/wrapped cell.
# Uses path keywords and font heuristics.
# ---------------------------------------------------------------------------
.is_title_like <- function(tg) {
  path <- tolower(tg$path)
  fs   <- if (!is.na(tg$fontsize)) tg$fontsize else 12

  # Explicit title/subtitle cells in nested gtables
  if (grepl("title", path) && !grepl("axis", path)) return(TRUE)

  # Heuristic: >=12pt text longer than 3 chars, not an axis label
  if (fs >= 12 && nchar(tg$label) > 3 && !grepl("axis", path)) return(TRUE)

  FALSE
}

# ---------------------------------------------------------------------------
# Helper: extract panel number from cell name ("tag-2" -> 2, "full-3" -> 3)
# ---------------------------------------------------------------------------
.cell_number <- function(name) {
  m <- regmatches(name, regexpr("[0-9]+$", name))
  if (length(m) == 1) as.integer(m) else NA_integer_
}

# ---------------------------------------------------------------------------
# Helper: issue constructor
# ---------------------------------------------------------------------------
.make_issue <- function(check, severity, message, elements = list(),
                         suggestion = "") {
  list(check = check, severity = severity, message = message,
       elements = elements, suggestion = suggestion)
}

# ===========================================================================
# Main validation function
# ===========================================================================
validate_composite_layout <- function(
  p,
  fig_width  = 16,
  fig_height = 18,
  # A title at npc y >= this threshold is considered "at the top" of its cell,
  # meaning it could collide with a tag in the cell above.
  top_threshold = 0.85,
  verbose    = TRUE
) {

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork package is required")

  errors   <- list()
  warnings <- list()
  info     <- list()

  # -----------------------------------------------------------------------
  # Phase A: Build gtable, compute cell positions in inches
  # -----------------------------------------------------------------------
  gt <- patchwork::patchworkGrob(p)

  # We need a graphics device to convert grid units to inches
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = fig_width, height = fig_height,
                 units = "in", res = 72)
  on.exit({ try(grDevices::dev.off(), silent = TRUE); unlink(tmp) })
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(fig_width, "in"),
    height = grid::unit(fig_height, "in")
  ))
  w_in <- grid::convertWidth(gt$widths, "inches", valueOnly = TRUE)
  h_in <- grid::convertHeight(gt$heights, "inches", valueOnly = TRUE)
  grid::popViewport()
  grDevices::dev.off()
  on.exit()  # clear the on.exit now that device is closed
  unlink(tmp)

  cum_w <- c(0, cumsum(w_in))  # left edges of columns (+ right edge of last)
  cum_h <- c(0, cumsum(h_in))  # top edges of rows (+ bottom of last)

  # -----------------------------------------------------------------------
  # Phase B: Catalog layout cells by type
  # -----------------------------------------------------------------------
  # tag_cells:  list of (panel_num, row_t, row_b, col_l, col_r, grob_idx)
  # panel_cells: same structure for full-N, panel-N, panel_area-N cells
  # title_cells: same for title-N cells
  tag_cells   <- list()
  panel_cells <- list()
  title_cells <- list()

  for (i in seq_len(nrow(gt$layout))) {
    nm  <- gt$layout$name[i]
    pn  <- .cell_number(nm)
    rec <- list(name = nm, panel_num = pn,
                t = gt$layout$t[i], b = gt$layout$b[i],
                l = gt$layout$l[i], r = gt$layout$r[i],
                grob_idx = i)

    if (grepl("^tag-", nm))                                tag_cells[[length(tag_cells) + 1]] <- rec
    else if (grepl("^(full|panel)-[0-9]+$", nm))           panel_cells[[length(panel_cells) + 1]] <- rec
    else if (grepl("^title-", nm) && !grepl("^axis", nm))  title_cells[[length(title_cells) + 1]] <- rec
  }

  # -----------------------------------------------------------------------
  # Phase C: Extract tag labels and their info
  # -----------------------------------------------------------------------
  tags <- list()  # list of (label, panel_num, cell_name, fontsize, fontface, bottom_y)
  for (tc in tag_cells) {
    grob <- gt$grobs[[tc$grob_idx]]
    texts <- .find_text_grobs(grob)
    for (tg in texts) {
      if (nchar(tg$label) == 1 && grepl("^[A-Z]$", tg$label)) {
        tags[[length(tags) + 1]] <- list(
          label     = tg$label,
          panel_num = tc$panel_num,
          cell      = tc$name,
          fontsize  = tg$fontsize,
          fontface  = tg$fontface,
          # Bottom edge of the tag cell in inches from top of figure
          bottom_y  = cum_h[tc$b + 1]
        )
      }
    }
  }

  n_tags <- length(tags)

  # -----------------------------------------------------------------------
  # Phase D: Find title-like text at the top of panel/wrapped cells
  # -----------------------------------------------------------------------
  top_titles <- list()  # list of (label, panel_num, cell_name, y_npc, top_y)
  for (pc in panel_cells) {
    grob <- gt$grobs[[pc$grob_idx]]

    # For wrapped elements containing a nested patchwork, the dangerous case
    # is a plot_annotation(title=...) which creates a cell named exactly
    # "title" (no numeric suffix) at the very top of the inner gtable.
    # Individual plot titles are named "title-N" and are safely nested
    # within their own plot's layout, away from the outer tag.
    if (inherits(grob, "gtable")) {
      for (gi in seq_len(nrow(grob$layout))) {
        cell_name <- grob$layout$name[gi]
        # Match "title" or "subtitle" exactly (patchwork annotation level),
        # but NOT "title-1", "title-2", etc. (individual plot level)
        if (cell_name %in% c("title", "subtitle")) {
          inner_grob <- grob$grobs[[gi]]
          inner_texts <- .find_text_grobs(inner_grob)
          for (tg in inner_texts) {
            if (nchar(tg$label) > 0) {
              top_titles[[length(top_titles) + 1]] <- list(
                label     = tg$label,
                panel_num = pc$panel_num,
                cell      = pc$name,
                inner_cell = cell_name,
                y_npc     = tg$y_npc,
                top_y     = cum_h[pc$t]
              )
            }
          }
        }
      }
    } else {
      # Non-gtable wrapped element: fall back to scanning all text
      texts <- .find_text_grobs(grob)
      for (tg in texts) {
        if (.is_title_like(tg) && !is.na(tg$y_npc) && tg$y_npc >= top_threshold) {
          top_titles[[length(top_titles) + 1]] <- list(
            label     = tg$label,
            panel_num = pc$panel_num,
            cell      = pc$name,
            inner_cell = "unknown",
            y_npc     = tg$y_npc,
            top_y     = cum_h[pc$t]
          )
        }
      }
    }
  }

  # Also count titles in their own layout cells (these are safe, just counted)
  n_safe_titles <- 0
  for (tc in title_cells) {
    grob <- gt$grobs[[tc$grob_idx]]
    texts <- .find_text_grobs(grob)
    n_safe_titles <- n_safe_titles + sum(vapply(texts,
      function(t) nchar(t$label) > 0, logical(1)))
  }

  n_titles <- length(top_titles) + n_safe_titles

  info[[length(info) + 1]] <- .make_issue(
    "composite_summary", "INFO",
    sprintf("Found %d tags, %d titles (%d in own cells, %d at top of panel cells)",
            n_tags, n_titles, n_safe_titles, length(top_titles)))

  # -----------------------------------------------------------------------
  # Check 1: Tag-title collision
  # A tag in tag-N can collide with a title at the top of the panel cell
  # below it (same column). The collision happens when the title is rendered
  # at npc y ≈ 1.0 in the panel cell, placing it right at the boundary
  # with the tag cell above.
  # -----------------------------------------------------------------------
  for (tag in tags) {
    for (title in top_titles) {
      # Must be in the same panel (same column group)
      if (!is.na(tag$panel_num) && !is.na(title$panel_num) &&
          tag$panel_num != title$panel_num) next

      # The tag's bottom edge should be close to the title's cell top edge
      gap <- abs(tag$bottom_y - title$top_y)
      gap_threshold <- fig_height * 0.03  # ~3% of figure height
      if (gap < gap_threshold) {  # likely adjacent cells
        errors[[length(errors) + 1]] <- .make_issue(
          "tag_title_overlap", "ERROR",
          sprintf("Panel tag '%s' likely overlaps title '%s' (title at y_npc=%.2f in %s, gap=%.2fin)",
                  tag$label, substr(title$label, 1, 40), title$y_npc, title$cell, gap),
          list(tag = tag$label, title = title$label,
               tag_cell = tag$cell, title_cell = title$cell),
          "Move the title from plot_annotation() into labs(title=...) on the inner plot, or add top margin"
        )
      }
    }
  }

  # -----------------------------------------------------------------------
  # Check 2: Tag-tag overlap
  # Two tags in the same row/column range would overlap.
  # -----------------------------------------------------------------------
  if (n_tags > 1) {
    for (i in seq_len(n_tags - 1)) {
      for (j in (i + 1):n_tags) {
        ti <- tags[[i]]
        tj <- tags[[j]]
        # Find their cells in tag_cells to compare row/col
        if (is.na(ti$panel_num) || is.na(tj$panel_num)) next
        idx_i <- which(vapply(tag_cells,
          function(tc) !is.na(tc$panel_num) && tc$panel_num == ti$panel_num, logical(1)))
        idx_j <- which(vapply(tag_cells,
          function(tc) !is.na(tc$panel_num) && tc$panel_num == tj$panel_num, logical(1)))
        if (length(idx_i) == 0 || length(idx_j) == 0) next
        tci <- tag_cells[[idx_i[1]]]
        tcj <- tag_cells[[idx_j[1]]]

        row_overlap <- tci$t <= tcj$b && tci$b >= tcj$t
        col_overlap <- tci$l <= tcj$r && tci$r >= tcj$l

        if (row_overlap && col_overlap) {
          errors[[length(errors) + 1]] <- .make_issue(
            "tag_tag_overlap", "ERROR",
            sprintf("Panel tags '%s' and '%s' occupy overlapping cells",
                    ti$label, tj$label),
            list(tag_1 = ti$label, tag_2 = tj$label),
            "Adjust layout so each tag has its own cell"
          )
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # Check 3: Missing or duplicate tags
  # -----------------------------------------------------------------------
  if (n_tags > 0) {
    tag_labels <- sort(vapply(tags, function(t) t$label, character(1)))

    # Duplicates
    dup <- tag_labels[duplicated(tag_labels)]
    if (length(dup) > 0) {
      errors[[length(errors) + 1]] <- .make_issue(
        "duplicate_tags", "ERROR",
        sprintf("Duplicate panel tags: %s", paste(unique(dup), collapse = ", ")),
        list(duplicates = unique(dup)),
        "Each panel should have a unique tag letter"
      )
    }

    # Gaps — only warn, since wrapped elements may contain their own tags
    unique_tags <- unique(tag_labels)
    max_pos <- max(match(unique_tags, LETTERS), na.rm = TRUE)
    if (!is.na(max_pos)) {
      expected <- LETTERS[seq_len(max_pos)]
      missing  <- setdiff(expected, unique_tags)
      if (length(missing) > 0) {
        warnings[[length(warnings) + 1]] <- .make_issue(
          "missing_tags", "WARNING",
          sprintf("Tag sequence {%s} has gaps: missing {%s}. Expected if tags are inside wrap_elements().",
                  paste(unique_tags, collapse = ","),
                  paste(missing, collapse = ",")),
          list(found = unique_tags, missing = missing),
          "Use manual tag_levels if some tags are inside wrapped panels"
        )
      }
    }
  }

  # -----------------------------------------------------------------------
  # Check 4: Tag font consistency
  # -----------------------------------------------------------------------
  if (n_tags > 1) {
    tag_sizes <- vapply(tags, function(t) {
      if (!is.na(t$fontsize)) t$fontsize else 0
    }, numeric(1))
    tag_sizes <- tag_sizes[tag_sizes > 0]

    if (length(unique(tag_sizes)) > 1) {
      warnings[[length(warnings) + 1]] <- .make_issue(
        "tag_font_inconsistency", "WARNING",
        sprintf("Tags use different font sizes: %s",
                paste(paste0(unique(tag_sizes), "pt"), collapse = ", ")),
        list(sizes = unique(tag_sizes)),
        "Use a consistent font size for all panel tags"
      )
    }
  }

  # -----------------------------------------------------------------------
  # Report
  # -----------------------------------------------------------------------
  result <- list(
    errors   = errors,
    warnings = warnings,
    info     = info,
    tags     = tags,
    top_titles = top_titles,
    summary  = list(
      n_tags     = n_tags,
      n_titles   = n_titles,
      n_errors   = length(errors),
      n_warnings = length(warnings),
      n_info     = length(info)
    )
  )

  if (verbose) {
    cat("\n=== Composite Figure Layout Validation ===\n")
    cat(sprintf("Tags: %d | Titles: %d (%d safe, %d near-top)\n\n",
                n_tags, n_titles, n_safe_titles, length(top_titles)))

    if (length(errors) > 0) {
      cat(sprintf("ERRORS (%d):\n", length(errors)))
      for (i in seq_along(errors)) {
        cat(sprintf("  [%d] %s: %s\n", i, errors[[i]]$check, errors[[i]]$message))
        if (nzchar(errors[[i]]$suggestion))
          cat(sprintf("      -> %s\n", errors[[i]]$suggestion))
      }
      cat("\n")
    }

    if (length(warnings) > 0) {
      cat(sprintf("WARNINGS (%d):\n", length(warnings)))
      for (i in seq_along(warnings)) {
        cat(sprintf("  [%d] %s: %s\n", i, warnings[[i]]$check, warnings[[i]]$message))
        if (nzchar(warnings[[i]]$suggestion))
          cat(sprintf("      -> %s\n", warnings[[i]]$suggestion))
      }
      cat("\n")
    }

    if (length(info) > 0) {
      cat(sprintf("INFO (%d):\n", length(info)))
      for (i in seq_along(info)) {
        cat(sprintf("  [%d] %s: %s\n", i, info[[i]]$check, info[[i]]$message))
      }
      cat("\n")
    }

    if (length(errors) == 0 && length(warnings) == 0) {
      cat("All composite checks passed.\n")
    } else if (length(errors) == 0) {
      cat(sprintf("No errors. %d warning(s) to review.\n", length(warnings)))
    } else {
      cat(sprintf("Result: %d error(s) must be fixed.\n", length(errors)))
    }
    cat("\n")
  }

  invisible(result)
}
