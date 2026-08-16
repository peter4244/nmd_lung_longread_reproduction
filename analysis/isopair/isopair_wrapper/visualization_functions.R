#!/usr/bin/env Rscript
#
# UCSC-Style Isoform Visualization Functions
#
# Publication-quality isoform pair diagrams in the style of the UCSC Genome
# Browser: thick rectangles for CDS, thin rectangles for UTR, directional
# arrows in introns.
#
# Requires: ggplot2

# ==============================================================================
# Data Building Functions
# ==============================================================================

#' Segment exons by CDS boundaries into UTR/CDS regions
#'
#' @param exon_starts Integer vector of exon start positions
#' @param exon_ends Integer vector of exon end positions
#' @param cds_start Integer; CDS start (genomic, min coordinate)
#' @param cds_stop Integer; CDS stop (genomic, max coordinate)
#' @param strand Character; "+" or "-"
#' @return Data frame with exon_start, exon_end, region columns
segment_exons_by_cds <- function(exon_starts, exon_ends, cds_start, cds_stop,
                                  strand) {
  if (strand == "+") {
    before_cds <- "5UTR"
    after_cds  <- "3UTR"
  } else {
    before_cds <- "3UTR"
    after_cds  <- "5UTR"
  }

  segments <- vector("list", length(exon_starts) * 3L)
  k <- 0L
  for (i in seq_along(exon_starts)) {
    es <- exon_starts[i]
    ee <- exon_ends[i]

    if (ee < cds_start) {
      k <- k + 1L
      segments[[k]] <- data.frame(exon_start = es, exon_end = ee,
                                  region = before_cds)
      next
    }
    if (es > cds_stop) {
      k <- k + 1L
      segments[[k]] <- data.frame(exon_start = es, exon_end = ee,
                                  region = after_cds)
      next
    }
    # Portion before CDS
    if (es < cds_start) {
      k <- k + 1L
      segments[[k]] <- data.frame(exon_start = es, exon_end = cds_start - 1L,
                                  region = before_cds)
    }
    # CDS portion
    os <- max(es, cds_start)
    oe <- min(ee, cds_stop)
    if (os <= oe) {
      k <- k + 1L
      segments[[k]] <- data.frame(exon_start = os, exon_end = oe,
                                  region = "CDS")
    }
    # Portion after CDS
    if (ee > cds_stop) {
      k <- k + 1L
      segments[[k]] <- data.frame(exon_start = cds_stop + 1L, exon_end = ee,
                                  region = after_cds)
    }
  }
  do.call(rbind, segments[seq_len(k)])
}


#' Build rectangle data for a single isoform track (UCSC style)
#'
#' @param exon_starts Integer vector of exon start positions
#' @param exon_ends Integer vector of exon end positions
#' @param y_center Numeric; y-axis center for this track
#' @param strand Character; "+" or "-"
#' @param cds_start Integer; CDS start (or NULL if non-coding)
#' @param cds_stop Integer; CDS stop (or NULL if non-coding)
#' @param cds_height Numeric; half-height of CDS rectangles
#' @param utr_height Numeric; half-height of UTR rectangles
#' @return Data frame with xmin, xmax, ymin, ymax, region
build_ucsc_track <- function(exon_starts, exon_ends, y_center, strand,
                              cds_start = NULL, cds_stop = NULL,
                              cds_height = 0.30, utr_height = 0.12) {

  has_cds <- !is.null(cds_start) && !is.null(cds_stop)

  if (has_cds) {
    segs <- segment_exons_by_cds(exon_starts, exon_ends, cds_start, cds_stop,
                                  strand)
  } else {
    segs <- data.frame(exon_start = exon_starts, exon_end = exon_ends,
                       region = "non_coding")
  }

  segs$xmin <- segs$exon_start
  segs$xmax <- segs$exon_end
  segs$ymin <- ifelse(segs$region == "CDS",
                      y_center - cds_height, y_center - utr_height)
  segs$ymax <- ifelse(segs$region == "CDS",
                      y_center + cds_height, y_center + utr_height)
  segs$y_center <- y_center
  segs
}


#' Build intron line + directional arrow data
#'
#' @param exon_starts Integer vector of exon start positions (sorted ascending)
#' @param exon_ends Integer vector of exon end positions (sorted ascending)
#' @param y_center Numeric; y-axis center for this track
#' @param strand Character; "+" or "-"
#' @param arrow_size Numeric; half-height of arrowheads in y-units
#' @return List with $lines (data.frame) and $arrows (data.frame)
build_intron_arrows <- function(exon_starts, exon_ends, y_center, strand,
                                 arrow_size = 0.06,
                                 chev_spacing = NULL,
                                 chev_arm_width = NULL) {
  ord <- order(exon_starts)
  es <- exon_starts[ord]
  ee <- exon_ends[ord]
  n <- length(es)
  if (n < 2L) return(list(lines = data.frame(), arrows = data.frame()))

  # Intron line segments
  lines <- data.frame(
    x = ee[-n], xend = es[-1L],
    y = y_center, yend = y_center
  )

  # Use figure-level chevron parameters if provided; otherwise compute locally
  if (is.null(chev_spacing) || is.null(chev_arm_width)) {
    total_range <- max(ee) - min(es)
    chev_spacing <- total_range * 0.06
    chev_arm_width <- chev_spacing * 0.4
  }

  # Minimum intron length: must fit at least one chevron (arm + margins)
  min_intron <- chev_arm_width * 3

  ax <- numeric(0); axend <- numeric(0)
  ay <- numeric(0); ayend <- numeric(0)

  for (i in seq_len(n - 1L)) {
    intron_start <- ee[i]
    intron_end   <- es[i + 1L]
    intron_len   <- intron_end - intron_start
    if (intron_len < min_intron) next

    # Number of chevrons: at least 1, spaced by chev_spacing
    n_chev <- max(1L, floor(intron_len / chev_spacing))
    step <- intron_len / (n_chev + 1)
    positions <- intron_start + step * seq_len(n_chev)

    for (cx in positions) {
      if (strand == "+") {
        ax    <- c(ax,    cx - chev_arm_width, cx - chev_arm_width)
        axend <- c(axend, cx,                  cx)
      } else {
        ax    <- c(ax,    cx + chev_arm_width, cx + chev_arm_width)
        axend <- c(axend, cx,                  cx)
      }
      ay    <- c(ay,    y_center + arrow_size, y_center - arrow_size)
      ayend <- c(ayend, y_center,              y_center)
    }
  }

  arrows <- data.frame(x = ax, xend = axend, y = ay, yend = ayend)
  list(lines = lines, arrows = arrows)
}


# ==============================================================================
# Main Plotting Function
# ==============================================================================

#' Plot an isoform pair in UCSC Genome Browser style
#'
#' Thick CDS rectangles, thin UTR rectangles, directional intron arrows.
#'
#' @param ref_starts,ref_ends Integer vectors; reference exon coordinates
#' @param comp_starts,comp_ends Integer vectors; comparator exon coordinates
#' @param strand Character; "+" or "-"
#' @param ref_cds Named list with cds_start, cds_stop (or NULL)
#' @param comp_cds Named list with cds_start, cds_stop (or NULL)
#' @param ref_label Character; label for reference track
#' @param comp_label Character; label for comparator track
#' @param events Data frame from detectEvents() (optional, for annotation brackets)
#' @param annotations List of annotation specs (see details)
#' @param title Character; plot title
#' @param subtitle Character; plot subtitle
#' @param colors Named character vector overriding default region colors
#' @return A ggplot object
#'
#' @details
#' Annotations are a list of lists, each with:
#' \describe{
#'   \item{type}{"ptc_marker", "frameshift_boundary", or "region_shade"}
#'   \item{x}{Genomic x position (for ptc_marker and frameshift_boundary)}
#'   \item{track}{"ref" or "comp" (for ptc_marker)}
#'   \item{xmin, xmax}{Genomic range (for region_shade)}
#'   \item{track}{"ref" or "comp" (for region_shade)}
#'   \item{label}{Optional text label}
#'   \item{color}{Optional color override}
#' }
plot_ucsc_pair <- function(ref_starts, ref_ends,
                            comp_starts, comp_ends,
                            strand,
                            ref_cds = NULL, comp_cds = NULL,
                            ref_label = "Reference",
                            comp_label = "Comparator",
                            events = NULL,
                            annotations = NULL,
                            title = NULL, subtitle = NULL,
                            colors = NULL,
                            xlim = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required")

  # Default colors
  default_colors <- c(
    "CDS"        = "#2C3E50",
    "5UTR"       = "#7F8C8D",
    "3UTR"       = "#7F8C8D",
    "non_coding" = "#BDC3C7"
  )
  if (!is.null(colors)) default_colors[names(colors)] <- colors

  ref_y  <- 2

  comp_y <- 1

  # Build track data
  ref_track <- build_ucsc_track(
    ref_starts, ref_ends, ref_y, strand,
    ref_cds$cds_start, ref_cds$cds_stop
  )
  comp_track <- build_ucsc_track(
    comp_starts, comp_ends, comp_y, strand,
    comp_cds$cds_start, comp_cds$cds_stop
  )
  all_rects <- rbind(
    cbind(ref_track, track = "ref"),
    cbind(comp_track, track = "comp")
  )

  # Build intron data — compute figure-level chevron parameters so all
  # chevrons in the figure have identical size and spacing
  all_coords <- c(ref_starts, ref_ends, comp_starts, comp_ends)
  fig_range <- max(all_coords) - min(all_coords)
  fig_chev_spacing <- fig_range * 0.06
  fig_chev_arm <- fig_chev_spacing * 0.4

  ref_introns  <- build_intron_arrows(ref_starts, ref_ends, ref_y, strand,
                                       chev_spacing = fig_chev_spacing,
                                       chev_arm_width = fig_chev_arm)
  comp_introns <- build_intron_arrows(comp_starts, comp_ends, comp_y, strand,
                                       chev_spacing = fig_chev_spacing,
                                       chev_arm_width = fig_chev_arm)
  all_lines  <- rbind(ref_introns$lines, comp_introns$lines)
  all_arrows <- rbind(ref_introns$arrows, comp_introns$arrows)

  # Genomic range
  x_min <- min(all_rects$xmin)
  x_max <- max(all_rects$xmax)
  x_range <- x_max - x_min
  x_buf <- max(x_range * 0.05, 50)

  # Build plot
  p <- ggplot2::ggplot()

  # Intron lines
  if (nrow(all_lines) > 0) {
    p <- p + ggplot2::geom_segment(
      data = all_lines,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                   y = .data$y, yend = .data$yend),
      color = "gray40", linewidth = 0.4
    )
  }

  # Directional arrows
  if (nrow(all_arrows) > 0) {
    p <- p + ggplot2::geom_segment(
      data = all_arrows,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                   y = .data$y, yend = .data$yend),
      color = "gray40", linewidth = 0.3
    )
  }

  # Exon rectangles
  p <- p + ggplot2::geom_rect(
    data = all_rects,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                 ymin = .data$ymin, ymax = .data$ymax,
                 fill = .data$region),
    color = "black", linewidth = 0.3
  ) +
    ggplot2::scale_fill_manual(
      values = default_colors,
      name = NULL,
      breaks = c("CDS", "5UTR", "3UTR"),
      labels = c("CDS", "5\u2032 UTR", "3\u2032 UTR")
    )

  # Event annotation brackets
  if (!is.null(events) && nrow(events) > 0) {
    y_base <- 0.4
    y_step <- 0.28
    events <- events[order(events$five_prime), ]
    for (i in seq_len(nrow(events))) {
      ev <- events[i, ]
      y_level <- y_base - (i - 1) * y_step
      dir_sym <- if (ev$direction == "LOSS") "\u2193" else "\u2191"
      ev_label <- sprintf("%s %s", ev$event_type, dir_sym)
      br_color <- if (ev$direction == "LOSS") "#C0392B" else "#2980B9"

      p <- p +
        ggplot2::geom_segment(
          ggplot2::aes(x = !!ev$five_prime, xend = !!ev$three_prime,
                       y = !!y_level, yend = !!y_level),
          color = br_color, linewidth = 1.0
        ) +
        ggplot2::geom_segment(
          ggplot2::aes(x = !!ev$five_prime, xend = !!ev$five_prime,
                       y = !!y_level, yend = !!(y_level + 0.10)),
          color = br_color, linewidth = 0.7
        ) +
        ggplot2::geom_segment(
          ggplot2::aes(x = !!ev$three_prime, xend = !!ev$three_prime,
                       y = !!y_level, yend = !!(y_level + 0.10)),
          color = br_color, linewidth = 0.7
        ) +
        ggplot2::annotate("text",
          x = (ev$five_prime + ev$three_prime) / 2, y = y_level - 0.10,
          label = ev_label, size = 2.5, color = br_color, fontface = "bold"
        )
    }
  }

  # Custom annotations
  if (!is.null(annotations)) {
    for (ann in annotations) {
      ann_y <- if (ann$track == "ref") ref_y else comp_y

      if (ann$type == "ptc_marker") {
        marker_color <- if (!is.null(ann$color)) ann$color else "#E74C3C"
        p <- p +
          ggplot2::annotate("point",
            x = ann$x, y = ann_y + 0.42,
            shape = 25, size = 3, fill = marker_color, color = marker_color
          )
        if (!is.null(ann$label)) {
          p <- p + ggplot2::annotate("text",
            x = ann$x, y = ann_y + 0.55,
            label = ann$label, size = 2.3, color = marker_color,
            fontface = "bold"
          )
        }

      } else if (ann$type == "frameshift_boundary") {
        fs_color <- if (!is.null(ann$color)) ann$color else "#8E44AD"
        p <- p +
          ggplot2::geom_vline(xintercept = ann$x,
                              linetype = "dashed", color = fs_color,
                              linewidth = 0.6)
        if (!is.null(ann$label)) {
          p <- p + ggplot2::annotate("text",
            x = ann$x, y = max(ref_y, comp_y) + 0.55,
            label = ann$label, size = 2.3, color = fs_color,
            fontface = "bold", hjust = -0.1
          )
        }

      } else if (ann$type == "region_shade") {
        shade_color <- if (!is.null(ann$color)) ann$color else "#F39C12"
        shade_alpha <- if (!is.null(ann$alpha)) ann$alpha else 0.25
        p <- p +
          ggplot2::annotate("rect",
            xmin = ann$xmin, xmax = ann$xmax,
            ymin = ann_y - 0.40, ymax = ann_y + 0.40,
            fill = shade_color, alpha = shade_alpha
          )
        if (!is.null(ann$label)) {
          p <- p + ggplot2::annotate("text",
            x = (ann$xmin + ann$xmax) / 2, y = ann_y - 0.50,
            label = ann$label, size = 2.3, color = shade_color,
            fontface = "italic"
          )
        }
      }
    }
  }

  # Track labels
  label_df <- data.frame(
    y = c(ref_y, comp_y),
    label = c(ref_label, comp_label)
  )
  y_bottom <- if (!is.null(events) && nrow(events) > 1) {
    0.4 - nrow(events) * 0.28 - 0.3
  } else {
    0.2
  }

  p <- p +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x_min - x_buf * 1.2, y = .data$y, label = .data$label),
      hjust = 1, size = 3.2, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL,
                                limits = c(y_bottom, ref_y + 0.7)) +
    ggplot2::coord_cartesian(
      xlim = if (!is.null(xlim)) xlim else c(x_min - x_buf * 2.5, x_max + x_buf),
      clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "none",
      axis.text.y      = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(size = 8),
      plot.margin      = ggplot2::margin(5, 10, 5, 80),
      plot.title       = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "gray30")
    )

  if (!is.null(title))    p <- p + ggplot2::ggtitle(title, subtitle = subtitle)
  p
}


# ==============================================================================
# Flowchart / Funnel Diagram
# ==============================================================================

#' Build a box-and-arrow funnel diagram from counts
#'
#' @param stages Data frame with columns: label, n, x, y, width, height, fill
#' @param arrows Data frame with columns: from_label, to_label
#' @param title Character; plot title
#' @return A ggplot object
plot_funnel <- function(stages, arrows, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required")

  # Compute box boundaries
  stages$xmin <- stages$x - stages$width / 2
  stages$xmax <- stages$x + stages$width / 2
  stages$ymin <- stages$y - stages$height / 2
  stages$ymax <- stages$y + stages$height / 2
  stages$display <- paste0(stages$label, "\n(n=", formatC(stages$n, big.mark = ","), ")")

  # Build arrow segments: connect bottom of source to top of target
  arrow_segs <- do.call(rbind, lapply(seq_len(nrow(arrows)), function(i) {
    from <- stages[stages$label == arrows$from_label[i], ]
    to   <- stages[stages$label == arrows$to_label[i], ]
    if (nrow(from) == 0 || nrow(to) == 0) return(NULL)
    data.frame(
      x = from$x, y = from$ymin,
      xend = to$x, yend = to$ymax
    )
  }))

  p <- ggplot2::ggplot()

  # Arrows
  if (!is.null(arrow_segs) && nrow(arrow_segs) > 0) {
    p <- p + ggplot2::geom_segment(
      data = arrow_segs,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                   y = .data$y, yend = .data$yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "inches"),
                             type = "closed"),
      color = "gray40", linewidth = 0.5
    )
  }

  # Boxes
  p <- p +
    ggplot2::geom_rect(
      data = stages,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax,
                   fill = .data$fill),
      color = "black", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = stages,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$display),
      size = 3.2, lineheight = 0.85
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(15, 20, 15, 20)
    )

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


# ==============================================================================
# Helper: Extract exons from structures object for a given isoform
# ==============================================================================

#' Get exon start/end vectors from structures data frame
#'
#' @param isoform_id Character; isoform to look up
#' @param structures Data frame with isoform_id, exon_starts (list), exon_ends (list)
#' @return List with $starts and $ends (integer vectors), or NULL
get_exon_coords <- function(isoform_id, structures) {
  row <- structures[structures$isoform_id == isoform_id, ]
  if (nrow(row) == 0) return(NULL)
  list(starts = row$exon_starts[[1]], ends = row$exon_ends[[1]])
}

#' Get CDS range from cds data frame
#'
#' @param isoform_id Character; isoform to look up
#' @param cds_df Data frame with isoform_id, cds_start, cds_stop
#' @return List with $cds_start and $cds_stop, or NULL
get_cds_range <- function(isoform_id, cds_df) {
  row <- cds_df[cds_df$isoform_id == isoform_id & cds_df$coding_status == "coding", ]
  if (nrow(row) == 0) return(NULL)
  list(cds_start = row$cds_start[1], cds_stop = row$cds_stop[1])
}
