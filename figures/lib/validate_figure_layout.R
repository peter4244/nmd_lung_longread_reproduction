#!/usr/bin/env Rscript
# validate_figure_layout.R
#
# Validates spatial layout of ggplot2 annotation-based figures (architecture
# diagrams, flowcharts, etc.).  Checks for text-arrow collisions, text
# overlap, crowding, centering, alignment, and padding.
#
# For patchwork composite figures, use validate_composite_layout.R instead.
#
# Usage:
#   source("figures/lib/validate_figure_layout.R")
#   source("make_architecture_figure.R")
#   result <- validate_figure_layout(p, fig_width = 16, fig_height = 24)
#
# All thresholds are in data-coordinate units and can be tuned via arguments.

# ---------------------------------------------------------------------------
# Helper: extract annotation layers from a ggplot object
# ---------------------------------------------------------------------------
.extract_layers <- function(p) {
  texts    <- list()
  segments <- list()
  rects    <- list()

  for (i in seq_along(p$layers)) {
    layer <- p$layers[[i]]
    geom  <- layer$geom
    d     <- layer$data
    ap    <- layer$aes_params  # aesthetic params set via annotate()

    if (inherits(geom, "GeomText") || inherits(geom, "GeomLabel")) {
      label <- if ("label" %in% names(d)) as.character(d$label) else
               if ("label" %in% names(ap)) as.character(ap$label) else NA_character_
      sz    <- if ("size" %in% names(ap)) ap$size else
               if ("size" %in% names(d)) d$size else 3.88  # ggplot default
      hj    <- if ("hjust" %in% names(ap)) ap$hjust else
               if ("hjust" %in% names(d)) d$hjust else 0.5
      vj    <- if ("vjust" %in% names(ap)) ap$vjust else
               if ("vjust" %in% names(d)) d$vjust else 0.5
      al    <- if ("alpha" %in% names(ap)) ap$alpha else
               if ("alpha" %in% names(d)) d$alpha else 1

      texts[[length(texts) + 1]] <- data.frame(
        idx   = i,
        x     = d$x,
        y     = d$y,
        label = label,
        size  = sz,
        hjust = hj,
        vjust = vj,
        alpha = al,
        stringsAsFactors = FALSE
      )
    } else if (inherits(geom, "GeomSegment")) {
      # Detect whether the layer was constructed with an arrow= spec — used to
      # skip arrow-only checks (e.g. minimum-length) for plain tick segments,
      # range bars, brackets, etc. ggplot2 stashes the arrow in slightly
      # different places across versions, so probe both:
      has_arrow <- isTRUE(inherits(layer$geom_params$arrow, "arrow")) ||
                   isTRUE(inherits(layer$arrow,             "arrow")) ||
                   isTRUE(inherits(ap$arrow,                "arrow"))
      segments[[length(segments) + 1]] <- data.frame(
        idx       = i,
        x         = d$x,
        xend      = d$xend,
        y         = d$y,
        yend      = d$yend,
        has_arrow = has_arrow,
        stringsAsFactors = FALSE
      )
    } else if (inherits(geom, "GeomRect")) {
      fl <- if ("fill" %in% names(ap)) ap$fill else
            if ("fill" %in% names(d)) d$fill else NA_character_
      al <- if ("alpha" %in% names(ap)) ap$alpha else
            if ("alpha" %in% names(d)) d$alpha else 1

      rects[[length(rects) + 1]] <- data.frame(
        idx   = i,
        xmin  = d$xmin,
        xmax  = d$xmax,
        ymin  = d$ymin,
        ymax  = d$ymax,
        fill  = fl,
        alpha = al,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    texts    = if (length(texts))    do.call(rbind, texts)    else data.frame(),
    segments = if (length(segments)) do.call(rbind, segments) else data.frame(),
    rects    = if (length(rects))    do.call(rbind, rects)    else data.frame()
  )
}

# ---------------------------------------------------------------------------
# Helper: compute text bounding boxes in data coordinates
#
# ggplot2 size is in points (1 pt = 0.353 mm).  Average character width
# is roughly 0.5 * size_mm; height is roughly 0.7 * size_mm (with leading).
# We convert physical mm to data-coordinate units using the figure dimensions
# and coordinate ranges.
# ---------------------------------------------------------------------------
.compute_text_bboxes <- function(texts, fig_width, fig_height,
                                  x_range, y_range, bbox_pad = 1.1) {
  if (nrow(texts) == 0) return(texts)

  # mm per data unit along each axis
  mm_per_unit_x <- (fig_width  * 25.4) / x_range
  mm_per_unit_y <- (fig_height * 25.4) / y_range

  texts$bb_xmin <- NA_real_
  texts$bb_xmax <- NA_real_
  texts$bb_ymin <- NA_real_
  texts$bb_ymax <- NA_real_

  for (i in seq_len(nrow(texts))) {
    lab  <- texts$label[i]
    sz   <- texts$size[i]
    hj   <- texts$hjust[i]
    vj   <- texts$vjust[i]

    # Size in mm (ggplot2 size unit = pt; 1 pt = 0.353 mm)
    sz_mm <- sz * 0.353

    # Handle multi-line labels
    lines   <- strsplit(lab, "\n", fixed = TRUE)[[1]]
    n_lines <- length(lines)
    max_nc  <- max(nchar(lines))

    # Approximate physical dimensions
    char_w_mm <- sz_mm * 0.55    # average character width
    char_h_mm <- sz_mm * 1.2     # line height with leading

    total_w_mm <- max_nc * char_w_mm
    total_h_mm <- n_lines * char_h_mm

    # Convert to data units with padding factor
    bbox_w <- (total_w_mm / mm_per_unit_x) * bbox_pad
    bbox_h <- (total_h_mm / mm_per_unit_y) * bbox_pad

    x <- texts$x[i]
    y <- texts$y[i]

    texts$bb_xmin[i] <- x - hj * bbox_w
    texts$bb_xmax[i] <- x + (1 - hj) * bbox_w
    texts$bb_ymin[i] <- y - vj * bbox_h
    texts$bb_ymax[i] <- y + (1 - vj) * bbox_h
  }

  texts
}

# ---------------------------------------------------------------------------
# Helper: Liang-Barsky line-rectangle intersection test
# Returns TRUE if segment (x1,y1)->(x2,y2) intersects rect [xmin,xmax,ymin,ymax]
# ---------------------------------------------------------------------------
.segment_intersects_rect <- function(x1, y1, x2, y2,
                                      xmin, xmax, ymin, ymax) {
  dx <- x2 - x1
  dy <- y2 - y1

  p <- c(-dx, dx, -dy, dy)
  q <- c(x1 - xmin, xmax - x1, y1 - ymin, ymax - y1)

  t0 <- 0
  t1 <- 1

  for (k in 1:4) {
    if (abs(p[k]) < 1e-12) {
      # Parallel to this edge
      if (q[k] < 0) return(FALSE)  # outside
    } else {
      r <- q[k] / p[k]
      if (p[k] < 0) {
        t0 <- max(t0, r)
      } else {
        t1 <- min(t1, r)
      }
      if (t0 > t1) return(FALSE)
    }
  }

  TRUE
}

# ---------------------------------------------------------------------------
# Helper: minimum distance from point to line segment
# ---------------------------------------------------------------------------
.point_to_segment_dist <- function(px, py, x1, y1, x2, y2) {
  dx <- x2 - x1
  dy <- y2 - y1
  len_sq <- dx^2 + dy^2

  if (len_sq < 1e-12) {
    return(sqrt((px - x1)^2 + (py - y1)^2))
  }

  t <- ((px - x1) * dx + (py - y1) * dy) / len_sq
  t <- max(0, min(1, t))

  proj_x <- x1 + t * dx
  proj_y <- y1 + t * dy

  sqrt((px - proj_x)^2 + (py - proj_y)^2)
}

# ---------------------------------------------------------------------------
# Helper: minimum distance from rectangle boundary to line segment
# Uses the 4 corners + 4 edge midpoints as sampling approximation
# ---------------------------------------------------------------------------
.rect_to_segment_dist <- function(rxmin, rxmax, rymin, rymax,
                                   sx1, sy1, sx2, sy2) {
  # Sample points on the rectangle perimeter
  pts_x <- c(rxmin, rxmax, rxmin, rxmax,
              (rxmin + rxmax) / 2, (rxmin + rxmax) / 2, rxmin, rxmax)
  pts_y <- c(rymin, rymin, rymax, rymax,
              rymin, rymax, (rymin + rymax) / 2, (rymin + rymax) / 2)

  dists <- vapply(seq_along(pts_x), function(k) {
    .point_to_segment_dist(pts_x[k], pts_y[k], sx1, sy1, sx2, sy2)
  }, numeric(1))

  # Also check segment endpoints to rect boundary
  for (ep in list(c(sx1, sy1), c(sx2, sy2))) {
    # Distance from endpoint to nearest rect edge
    dx <- max(rxmin - ep[1], 0, ep[1] - rxmax)
    dy <- max(rymin - ep[2], 0, ep[2] - rymax)
    dists <- c(dists, sqrt(dx^2 + dy^2))
  }

  min(dists)
}

# ---------------------------------------------------------------------------
# Helper: do two axis-aligned rects overlap?
# ---------------------------------------------------------------------------
.rects_overlap <- function(xmin1, xmax1, ymin1, ymax1,
                            xmin2, xmax2, ymin2, ymax2) {
  (xmin1 < xmax2) && (xmax1 > xmin2) && (ymin1 < ymax2) && (ymax1 > ymin2)
}

# ---------------------------------------------------------------------------
# Helper: is this a (text, segment) intentional label pair?
#
# A text positioned cleanly above OR below a segment whose x-range it
# overlaps is *labeling* that segment, not being crowded by it. Returns
# TRUE when both conditions hold:
#   1. text bbox x-range overlaps segment x-range
#   2. text bbox is entirely above or entirely below the segment (no y
#      overlap — the text doesn't pierce the segment's vertical extent)
#
# Real crowding (label adjacent to an UNRELATED segment) violates one of
# these: either the segment is off to the side (no x overlap) or the text
# overlaps the segment in y. Those cases still trip the proximity warning.
# ---------------------------------------------------------------------------
.is_label_pair <- function(tx_xmin, tx_xmax, tx_ymin, tx_ymax,
                            sx1, sy1, sx2, sy2) {
  seg_xlo <- min(sx1, sx2); seg_xhi <- max(sx1, sx2)
  seg_ylo <- min(sy1, sy2); seg_yhi <- max(sy1, sy2)
  x_overlap <- !(tx_xmax < seg_xlo || tx_xmin > seg_xhi)
  y_above   <- tx_ymin >= seg_yhi
  y_below   <- tx_ymax <= seg_ylo
  x_overlap && (y_above || y_below)
}

# ---------------------------------------------------------------------------
# Helper: create an issue record
# ---------------------------------------------------------------------------
.make_issue <- function(check, severity, message, elements = list(),
                         suggestion = "") {
  list(check = check, severity = severity, message = message,
       elements = elements, suggestion = suggestion)
}

# ---------------------------------------------------------------------------
# Helper: print formatted report
# ---------------------------------------------------------------------------
.print_report <- function(result) {
  s <- result$summary

  cat("\n=== Figure Layout Validation ===\n")
  cat(sprintf("Extracted: %d texts, %d segments, %d rects\n\n",
              s$n_texts, s$n_segments, s$n_rects))

  if (s$n_errors > 0) {
    cat(sprintf("ERRORS (%d):\n", s$n_errors))
    for (i in seq_along(result$errors)) {
      iss <- result$errors[[i]]
      cat(sprintf("  [%d] %s: %s\n", i, iss$check, iss$message))
      if (nzchar(iss$suggestion))
        cat(sprintf("      -> %s\n", iss$suggestion))
    }
    cat("\n")
  }

  if (s$n_warnings > 0) {
    cat(sprintf("WARNINGS (%d):\n", s$n_warnings))
    for (i in seq_along(result$warnings)) {
      iss <- result$warnings[[i]]
      cat(sprintf("  [%d] %s: %s\n", i, iss$check, iss$message))
      if (nzchar(iss$suggestion))
        cat(sprintf("      -> %s\n", iss$suggestion))
    }
    cat("\n")
  }

  if (s$n_info > 0) {
    cat(sprintf("INFO (%d):\n", s$n_info))
    for (i in seq_along(result$info)) {
      iss <- result$info[[i]]
      cat(sprintf("  [%d] %s: %s\n", i, iss$check, iss$message))
    }
    cat("\n")
  }

  if (s$n_errors == 0 && s$n_warnings == 0) {
    cat("All checks passed.\n")
  } else if (s$n_errors == 0) {
    cat(sprintf("No errors. %d warning(s) to review.\n", s$n_warnings))
  } else {
    cat(sprintf("Result: %d error(s) must be fixed.\n", s$n_errors))
  }
  cat("\n")
}

# ===========================================================================
# Main validation function
# ===========================================================================
validate_figure_layout <- function(
  p,
  fig_width           = 16,
  fig_height          = 24,
  # Bounding-box tuning
  bbox_pad            = 1.1,    # multiplicative padding on estimated text bbox
  # Check thresholds (all in data-coordinate units)
  collision_buffer     = 0.05,
  crowding_distance    = 0.4,
  centering_tolerance  = 0.3,
  alignment_y_tol      = 0.2,
  alignment_x_tol      = 0.5,
  min_arrow_length     = 0.3,
  rect_padding_min     = 0.15,
  verbose              = TRUE
) {

  errors   <- list()
  warnings <- list()
  info     <- list()

  # -----------------------------------------------------------------------
  # Phase A: Extract layers and coordinate system
  # -----------------------------------------------------------------------
  layers <- .extract_layers(p)
  texts  <- layers$texts
  segs   <- layers$segments
  rects  <- layers$rects

  # Extract coordinate limits from the plot
  coord   <- p$coordinates
  x_range <- if (!is.null(coord$limits$x)) diff(coord$limits$x) else 20
  y_range <- if (!is.null(coord$limits$y)) diff(coord$limits$y) else 20

  # -----------------------------------------------------------------------
  # Phase B: Compute text bounding boxes
  # -----------------------------------------------------------------------
  texts <- .compute_text_bboxes(texts, fig_width, fig_height,
                                 x_range, y_range, bbox_pad)

  # -----------------------------------------------------------------------
  # Phase C: Run checks
  # -----------------------------------------------------------------------

  # --- Check 1: Text-segment collision (ERROR) ---------------------------
  if (nrow(texts) > 0 && nrow(segs) > 0) {
    for (ti in seq_len(nrow(texts))) {
      # Skip near-invisible text
      if (!is.na(texts$alpha[ti]) && texts$alpha[ti] < 0.1) next

      bb <- list(
        xmin = texts$bb_xmin[ti] - collision_buffer,
        xmax = texts$bb_xmax[ti] + collision_buffer,
        ymin = texts$bb_ymin[ti] - collision_buffer,
        ymax = texts$bb_ymax[ti] + collision_buffer
      )

      for (si in seq_len(nrow(segs))) {
        hit <- .segment_intersects_rect(
          segs$x[si], segs$y[si], segs$xend[si], segs$yend[si],
          bb$xmin, bb$xmax, bb$ymin, bb$ymax
        )
        if (hit) {
          lab_short <- substr(texts$label[ti], 1, 40)
          errors[[length(errors) + 1]] <- .make_issue(
            "text_segment_collision", "ERROR",
            sprintf("Text '%s' at (%.1f, %.1f) intersects segment (%.1f,%.1f)->(%.1f,%.1f)",
                    lab_short, texts$x[ti], texts$y[ti],
                    segs$x[si], segs$y[si], segs$xend[si], segs$yend[si]),
            list(text_idx = texts$idx[ti], seg_idx = segs$idx[si]),
            "Move text away from the arrow path or reroute the segment"
          )
        }
      }
    }
  }

  # --- Check 2: Text-text overlap (ERROR) --------------------------------
  if (nrow(texts) > 1) {
    pairs <- combn(nrow(texts), 2)
    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]
      j <- pairs[2, k]
      # Skip near-invisible text
      if ((!is.na(texts$alpha[i]) && texts$alpha[i] < 0.1) ||
          (!is.na(texts$alpha[j]) && texts$alpha[j] < 0.1)) next

      if (.rects_overlap(texts$bb_xmin[i], texts$bb_xmax[i],
                         texts$bb_ymin[i], texts$bb_ymax[i],
                         texts$bb_xmin[j], texts$bb_xmax[j],
                         texts$bb_ymin[j], texts$bb_ymax[j])) {
        lab_i <- substr(texts$label[i], 1, 30)
        lab_j <- substr(texts$label[j], 1, 30)
        errors[[length(errors) + 1]] <- .make_issue(
          "text_text_overlap", "ERROR",
          sprintf("Text '%s' at (%.1f,%.1f) overlaps '%s' at (%.1f,%.1f)",
                  lab_i, texts$x[i], texts$y[i],
                  lab_j, texts$x[j], texts$y[j]),
          list(text_idx_1 = texts$idx[i], text_idx_2 = texts$idx[j]),
          "Increase vertical or horizontal spacing between these labels"
        )
      }
    }
  }

  # --- Check 3: Text-segment proximity (WARNING) -------------------------
  # Skip "intentional label-pair" proximity: a text positioned cleanly above
  # or below a segment whose x-range it overlaps is labeling that segment,
  # not being crowded by it. A real crowding case (label too close to an
  # UNRELATED segment) violates one of those conditions — the segment is
  # off to the side (x ranges don't overlap), or the text overlaps the
  # segment in y (vertical collision, not labeling).
  if (nrow(texts) > 0 && nrow(segs) > 0) {
    for (ti in seq_len(nrow(texts))) {
      if (!is.na(texts$alpha[ti]) && texts$alpha[ti] < 0.1) next

      for (si in seq_len(nrow(segs))) {
        dist <- .rect_to_segment_dist(
          texts$bb_xmin[ti], texts$bb_xmax[ti],
          texts$bb_ymin[ti], texts$bb_ymax[ti],
          segs$x[si], segs$y[si], segs$xend[si], segs$yend[si]
        )
        if (dist > 0 && dist < crowding_distance) {
          if (.is_label_pair(texts$bb_xmin[ti], texts$bb_xmax[ti],
                             texts$bb_ymin[ti], texts$bb_ymax[ti],
                             segs$x[si], segs$y[si], segs$xend[si], segs$yend[si])) {
            next  # intentional label-segment pair; not a crowding violation
          }
          lab_short <- substr(texts$label[ti], 1, 30)
          warnings[[length(warnings) + 1]] <- .make_issue(
            "text_segment_proximity", "WARNING",
            sprintf("Text '%s' at (%.1f,%.1f) is %.2f units from segment (%.1f,%.1f)->(%.1f,%.1f)",
                    lab_short, texts$x[ti], texts$y[ti], dist,
                    segs$x[si], segs$y[si], segs$xend[si], segs$yend[si]),
            list(text_idx = texts$idx[ti], seg_idx = segs$idx[si],
                 distance = dist),
            "Consider increasing spacing to avoid visual crowding"
          )
        }
      }
    }
  }

  # --- Check 4: Text-rect centering (WARNING) ----------------------------
  # For each text, find the smallest enclosing rect (not every container)
  if (nrow(texts) > 0 && nrow(rects) > 0) {
    for (ti in seq_len(nrow(texts))) {
      if (!is.na(texts$alpha[ti]) && texts$alpha[ti] < 0.1) next

      tx <- texts$x[ti]
      ty <- texts$y[ti]

      # Find all rects that contain or are near the text
      best_ri   <- NA
      best_area <- Inf

      for (ri in seq_len(nrow(rects))) {
        # Skip border-only rects (fill = NA) — these are containers/outlines,
        # not boxes that text should be centered within
        if (is.na(rects$fill[ri])) next

        inside_x <- tx >= rects$xmin[ri] && tx <= rects$xmax[ri]
        near_y   <- ty >= (rects$ymin[ri] - 1.0) && ty <= (rects$ymax[ri] + 1.0)

        if (inside_x && near_y) {
          area <- (rects$xmax[ri] - rects$xmin[ri]) *
                  (rects$ymax[ri] - rects$ymin[ri])
          if (area < best_area) {
            best_area <- area
            best_ri   <- ri
          }
        }
      }

      if (!is.na(best_ri)) {
        ri <- best_ri
        rx_center <- (rects$xmin[ri] + rects$xmax[ri]) / 2
        offset <- abs(tx - rx_center)
        if (offset > centering_tolerance) {
          lab_short <- substr(texts$label[ti], 1, 30)
          warnings[[length(warnings) + 1]] <- .make_issue(
            "text_rect_centering", "WARNING",
            sprintf("Text '%s' at x=%.1f is %.2f units off-center from rect [%.1f, %.1f] (center=%.1f)",
                    lab_short, tx, offset, rects$xmin[ri], rects$xmax[ri], rx_center),
            list(text_idx = texts$idx[ti], rect_idx = rects$idx[ri],
                 offset = offset),
            sprintf("Adjust text x to %.1f", rx_center)
          )
        }
      }
    }
  }

  # --- Check 5: Horizontal alignment (INFO) ------------------------------
  if (nrow(texts) > 2) {
    # Group texts by approximate y-coordinate
    ys <- texts$y
    ord <- order(ys)
    groups <- list()
    current_group <- c(ord[1])

    for (k in 2:length(ord)) {
      if (abs(ys[ord[k]] - ys[ord[k - 1]]) <= alignment_y_tol) {
        current_group <- c(current_group, ord[k])
      } else {
        if (length(current_group) >= 3) {
          groups[[length(groups) + 1]] <- current_group
        }
        current_group <- c(ord[k])
      }
    }
    if (length(current_group) >= 3) {
      groups[[length(groups) + 1]] <- current_group
    }

    # Within each group, check for near-but-not-matching y-values
    for (g in groups) {
      group_ys <- ys[g]
      mean_y   <- mean(group_ys)
      deviants <- which(abs(group_ys - mean_y) > alignment_y_tol)

      for (d in deviants) {
        ti <- g[d]
        lab_short <- substr(texts$label[ti], 1, 30)
        info[[length(info) + 1]] <- .make_issue(
          "horizontal_alignment", "INFO",
          sprintf("Text '%s' at y=%.2f deviates from group mean y=%.2f (group of %d texts)",
                  lab_short, texts$y[ti], mean_y, length(g)),
          list(text_idx = texts$idx[ti]),
          sprintf("Consider aligning y to %.2f", mean_y)
        )
      }
    }
  }

  # --- Check 6: Arrow minimum length (WARNING) ---------------------------
  # Only applied to segments constructed with an explicit arrow= spec.
  # Plain tick marks, bracket caps, range lines, etc. are not "arrows" and
  # are allowed to be short.
  if (nrow(segs) > 0) {
    arrow_segs <- if ("has_arrow" %in% names(segs)) segs[segs$has_arrow, , drop = FALSE]
                  else segs
    for (si in seq_len(nrow(arrow_segs))) {
      len <- sqrt((arrow_segs$xend[si] - arrow_segs$x[si])^2 +
                  (arrow_segs$yend[si] - arrow_segs$y[si])^2)
      if (len < min_arrow_length) {
        warnings[[length(warnings) + 1]] <- .make_issue(
          "arrow_min_length", "WARNING",
          sprintf("Arrow segment (%.1f,%.1f)->(%.1f,%.1f) length=%.2f < minimum %.2f",
                  arrow_segs$x[si], arrow_segs$y[si],
                  arrow_segs$xend[si], arrow_segs$yend[si],
                  len, min_arrow_length),
          list(seg_idx = arrow_segs$idx[si], length = len),
          "Increase arrow length or convert to a plain tick (drop the arrow= spec)"
        )
      }
    }
  }

  # --- Check 7: Text-rect padding (WARNING/ERROR) ------------------------
  if (nrow(texts) > 0 && nrow(rects) > 0) {
    for (ti in seq_len(nrow(texts))) {
      if (!is.na(texts$alpha[ti]) && texts$alpha[ti] < 0.1) next

      for (ri in seq_len(nrow(rects))) {
        # Skip border-only rects (fill = NA)
        if (is.na(rects$fill[ri])) next

        # Is text bbox inside this rect?
        inside <- (texts$bb_xmin[ti] >= rects$xmin[ri] - 0.01) &&
                  (texts$bb_xmax[ti] <= rects$xmax[ri] + 0.01) &&
                  (texts$bb_ymin[ti] >= rects$ymin[ri] - 0.01) &&
                  (texts$bb_ymax[ti] <= rects$ymax[ri] + 0.01)

        if (!inside) {
          # Check for overflow: text center is inside but bbox extends outside
          center_inside <- (texts$x[ti] >= rects$xmin[ri]) &&
                           (texts$x[ti] <= rects$xmax[ri]) &&
                           (texts$y[ti] >= rects$ymin[ri]) &&
                           (texts$y[ti] <= rects$ymax[ri])

          if (center_inside) {
            overflows <- c()
            if (texts$bb_xmin[ti] < rects$xmin[ri])
              overflows <- c(overflows, "left")
            if (texts$bb_xmax[ti] > rects$xmax[ri])
              overflows <- c(overflows, "right")
            if (texts$bb_ymin[ti] < rects$ymin[ri])
              overflows <- c(overflows, "bottom")
            if (texts$bb_ymax[ti] > rects$ymax[ri])
              overflows <- c(overflows, "top")

            if (length(overflows) > 0) {
              lab_short <- substr(texts$label[ti], 1, 30)
              errors[[length(errors) + 1]] <- .make_issue(
                "text_rect_overflow", "ERROR",
                sprintf("Text '%s' overflows rect [%.1f,%.1f]x[%.1f,%.1f] on: %s",
                        lab_short, rects$xmin[ri], rects$xmax[ri],
                        rects$ymin[ri], rects$ymax[ri],
                        paste(overflows, collapse = ", ")),
                list(text_idx = texts$idx[ti], rect_idx = rects$idx[ri]),
                "Reduce text size or widen the rect"
              )
            }
          }
        } else {
          # Text is inside rect — check padding
          left_pad   <- texts$bb_xmin[ti] - rects$xmin[ri]
          right_pad  <- rects$xmax[ri] - texts$bb_xmax[ti]
          top_pad    <- rects$ymax[ri] - texts$bb_ymax[ti]
          bottom_pad <- texts$bb_ymin[ti] - rects$ymin[ri]
          min_pad    <- min(left_pad, right_pad, top_pad, bottom_pad)

          if (min_pad < rect_padding_min) {
            lab_short <- substr(texts$label[ti], 1, 30)
            tight_side <- c("left", "right", "top", "bottom")[
              which.min(c(left_pad, right_pad, top_pad, bottom_pad))
            ]
            warnings[[length(warnings) + 1]] <- .make_issue(
              "text_rect_padding", "WARNING",
              sprintf("Text '%s' has only %.2f units padding on %s edge of rect [%.1f,%.1f]x[%.1f,%.1f]",
                      lab_short, min_pad, tight_side,
                      rects$xmin[ri], rects$xmax[ri],
                      rects$ymin[ri], rects$ymax[ri]),
              list(text_idx = texts$idx[ti], rect_idx = rects$idx[ri],
                   min_padding = min_pad, tight_side = tight_side),
              "Increase rect size or reduce text size"
            )
          }
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # Phase D: Assemble and report
  # -----------------------------------------------------------------------
  result <- list(
    errors   = errors,
    warnings = warnings,
    info     = info,
    summary  = list(
      n_texts    = nrow(texts),
      n_segments = nrow(segs),
      n_rects    = nrow(rects),
      n_errors   = length(errors),
      n_warnings = length(warnings),
      n_info     = length(info)
    )
  )

  if (verbose) .print_report(result)

  invisible(result)
}
