# ==============================================================================
# analysis_functions.R
#
# Shared analysis functions for the NMD isoform report.
# Sourced by both 05_final_report_mashr.Rmd and precompute scripts.
# ==============================================================================

#' Attribute PTC-causing splice events to individual pairs
#'
#' For each PTC+ pair, identifies the specific splice event responsible for the
#' premature stop codon using a three-mechanism framework:
#'   1. Frameshift: analyzeFrameWalk() event shifts reading frame → new stop
#'   2. In-frame stop: splice event introduces in-frame stop codon (event
#'      coordinates contain the PTC position)
#'   3. 3'UTR splice: same stop codon, but downstream EJC repositioned by
#'      a splicing event (handled separately, not by this function)
#'
#' Split-codon handling: a stop codon can span a splice junction (2 nt from
#' one exon, 1 from the next). The function extends coordinate checks by ±2 bp
#' to catch these cases.
#'
#' @param pairs Data frame with columns: comparator_isoform_id,
#'   reference_isoform_id. Each row is one PTC+ pair.
#' @param fw_events Data frame from analyzeFrameWalk()$events with columns:
#'   reference_isoform_id, comparator_isoform_id, event_type, genomic_start,
#'   genomic_end, is_frameshift.
#' @param profiles Data frame from buildProfiles() with detailed_events list
#'   column and comparator_isoform_id column.
#' @param ptc_genomic_pos Named numeric vector: comparator_isoform_id →
#'   strand-aware genomic position of the premature stop codon.
#' @param atg_genomic_pos Named numeric vector: comparator_isoform_id →
#'   genomic position of the CDS start ATG (for region filtering).
#'   If NULL, no region filtering is applied (all frame walk events considered).
#' @param strand_vec Named character vector: comparator_isoform_id → strand.
#' @param is_frameshift_vec Named logical vector: comparator_isoform_id →
#'   whether the pair is in a frameshift frame_category. If NULL, determined
#'   from fw_events (any frameshift event for this pair → TRUE).
#'
#' @return Data frame with columns: comparator_id, mechanism, ptc_causing_event,
#'   attribution ("direct" or "unresolved").
attribute_ptc_events <- function(pairs,
                                  fw_events,
                                  profiles,
                                  ptc_genomic_pos,
                                  atg_genomic_pos = NULL,
                                  strand_vec,
                                  is_frameshift_vec = NULL) {

  # Key on ref::comp to handle cases where the same comparator appears
  # in multiple pairs with different references (and different events)
  prof_keys <- paste0(profiles$reference_isoform_id, "::",
                       profiles$comparator_isoform_id)
  prof_idx_lookup <- setNames(seq_len(nrow(profiles)), prof_keys)

  result_rows <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    comp_id <- pairs$comparator_isoform_id[i]
    ref_id  <- pairs$reference_isoform_id[i]
    stop_g  <- ptc_genomic_pos[comp_id]
    strand  <- strand_vec[comp_id]

    if (is.na(stop_g) || is.na(strand)) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "unresolved",
        ptc_causing_event = "no_stop_pos", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    # Determine if this pair is in a frameshift category
    if (!is.null(is_frameshift_vec)) {
      is_fs <- isTRUE(is_frameshift_vec[comp_id])
    } else {
      pair_fw <- fw_events[fw_events$reference_isoform_id == ref_id &
                             fw_events$comparator_isoform_id == comp_id, ]
      is_fs <- any(pair_fw$is_frameshift, na.rm = TRUE)
    }

    # --- Mechanism 1: Frameshift ---
    if (is_fs) {
      pair_fw <- fw_events[fw_events$reference_isoform_id == ref_id &
                             fw_events$comparator_isoform_id == comp_id, ]

      if (nrow(pair_fw) > 0) {
        # If ATG position provided, filter to events overlapping ATG-to-stop
        # region (+3 bp buffer for split-codon splice junctions)
        if (!is.null(atg_genomic_pos)) {
          atg_g <- atg_genomic_pos[comp_id]
          if (!is.na(atg_g)) {
            region_lo <- min(atg_g, stop_g) - 3L
            region_hi <- max(atg_g, stop_g) + 3L
            pair_fw <- pair_fw[pair_fw$genomic_start <= region_hi &
                                 pair_fw$genomic_end >= region_lo, ]
          }
        }

        fs_events <- pair_fw[pair_fw$is_frameshift == TRUE, ]
        if (nrow(fs_events) > 0) {
          # First frameshift event in transcription direction
          if (!is.na(strand) && strand == "+") {
            first_fs <- fs_events[which.min(fs_events$genomic_start), ]
          } else {
            first_fs <- fs_events[which.max(fs_events$genomic_end), ]
          }
          result_rows[[i]] <- data.frame(
            comparator_isoform_id = comp_id, mechanism = "Frameshift",
            ptc_causing_event = first_fs$event_type[1],
            attribution = "direct", stringsAsFactors = FALSE)
          next
        }
      }

      # Frameshift category but no frameshift event found
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "Frameshift",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    # --- Mechanism 2: In-frame stop ---
    pair_key <- paste0(ref_id, "::", comp_id)
    prof_idx <- prof_idx_lookup[pair_key]
    if (is.na(prof_idx)) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    de <- profiles$detailed_events[[prof_idx]]
    if (is.null(de) || nrow(de) == 0) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    # Check which event's coordinates span the PTC position
    # ±2 bp buffer for split-codon splice junctions
    ev_min <- pmin(de$five_prime, de$three_prime)
    ev_max <- pmax(de$five_prime, de$three_prime)
    containing <- de[stop_g >= (ev_min - 2L) & stop_g <= (ev_max + 2L), ]

    if (nrow(containing) > 0) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = containing$event_type[1],
        attribution = "direct", stringsAsFactors = FALSE)
    } else {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
    }
  }

  do.call(rbind, result_rows)
}


#' Attribute 3'UTR splice PTC-causing events
#'
#' For same-stop PTC+ pairs (or same/longer-ORF pairs with downstream EJCs),
#' the PTC arises from a splice event downstream of the stop codon that
#' repositions an EJC. This function finds the nearest downstream event.
#'
#' @param pairs Data frame with comparator_isoform_id, reference_isoform_id
#' @param profiles Data frame with detailed_events list-column
#' @param stop_genomic_pos Named numeric vector: comp_id → stop genomic position
#' @param strand_vec Named character vector: comp_id → strand
#' @return Data frame with columns: comparator_id, mechanism, ptc_causing_event,
#'   attribution
attribute_3utr_splice <- function(pairs, profiles, stop_genomic_pos, strand_vec) {
  # Key on ref::comp to handle duplicate comparators
  prof_keys <- paste0(profiles$reference_isoform_id, "::",
                       profiles$comparator_isoform_id)
  prof_idx_lookup <- setNames(seq_len(nrow(profiles)), prof_keys)
  result_rows <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    comp_id <- pairs$comparator_isoform_id[i]
    ref_id  <- pairs$reference_isoform_id[i]
    stop_g <- stop_genomic_pos[comp_id]
    strand <- strand_vec[comp_id]

    if (is.na(stop_g) || is.na(strand)) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    pair_key <- paste0(ref_id, "::", comp_id)
    prof_idx <- prof_idx_lookup[pair_key]
    if (is.na(prof_idx)) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    de <- profiles$detailed_events[[prof_idx]]
    if (is.null(de) || nrow(de) == 0) {
      result_rows[[i]] <- data.frame(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved",
        stringsAsFactors = FALSE)
      next
    }

    # Find events downstream of stop in transcription direction
    # Use minimum absolute distance from either event boundary to the stop
    # to handle events that straddle the stop codon position
    if (strand == "+") {
      utr3_ev <- de[de$five_prime >= stop_g | de$three_prime >= stop_g, ]
      if (nrow(utr3_ev) == 0) {
        result_rows[[i]] <- data.frame(
          comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
          ptc_causing_event = "Unresolved", attribution = "unresolved",
          stringsAsFactors = FALSE)
        next
      }
      dist_to_stop <- pmin(abs(utr3_ev$five_prime - stop_g),
                            abs(utr3_ev$three_prime - stop_g))
      nearest_idx <- which.min(dist_to_stop)
    } else {
      utr3_ev <- de[de$five_prime <= stop_g | de$three_prime <= stop_g, ]
      if (nrow(utr3_ev) == 0) {
        result_rows[[i]] <- data.frame(
          comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
          ptc_causing_event = "Unresolved", attribution = "unresolved",
          stringsAsFactors = FALSE)
        next
      }
      dist_to_stop <- pmin(abs(utr3_ev$five_prime - stop_g),
                            abs(utr3_ev$three_prime - stop_g))
      nearest_idx <- which.min(dist_to_stop)
    }

    result_rows[[i]] <- data.frame(
      comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
      ptc_causing_event = utr3_ev$event_type[nearest_idx],
      attribution = "direct", stringsAsFactors = FALSE)
  }

  do.call(rbind, result_rows)
}


# ==============================================================================
# Display helpers
# ==============================================================================

#' Shorten splice event type labels for display
#'
#' @param event_types Character vector of event type names
#' @param collapse_rare Minimum count to keep; types below this become "Other"
#' @return Character vector with shortened labels
shorten_event_labels <- function(event_types, collapse_rare = 0) {
  labels <- event_types
  labels[labels == "Missing_Internal"] <- "Missing Int."
  labels[labels == "Partial_IR_5"] <- "Partial IR 5'"
  labels[labels == "Partial_IR_3"] <- "Partial IR 3'"
  labels[labels == "Alt_TSS"] <- "Alt TSS"
  labels[labels == "Alt_TES"] <- "Alt TES"
  labels[labels == "IR_diff"] <- "IR diff"
  labels[labels == "IR_diff_3"] <- "IR diff 3'"
  labels[labels == "IR_diff_5"] <- "IR diff 5'"
  if (collapse_rare > 0) {
    tab <- table(labels)
    rare <- names(tab[tab < collapse_rare])
    labels[labels %in% rare] <- "Other"
  }
  labels
}


# ==============================================================================
# CDS / coordinate helpers
# ==============================================================================

#' Get strand-aware stop codon genomic position
#'
#' On + strand, the stop codon is at cds_stop (larger coordinate).
#' On - strand, the stop codon is at cds_start (smaller coordinate).
#'
#' @param isoform_ids Character vector of isoform IDs
#' @param cds Data frame with columns: isoform_id, cds_start, cds_stop, strand
#' @return Named numeric vector of stop positions
get_strand_aware_stop <- function(isoform_ids, cds) {
  idx <- match(isoform_ids, cds$isoform_id)
  stops <- ifelse(cds$strand[idx] == "-", cds$cds_start[idx], cds$cds_stop[idx])
  names(stops) <- isoform_ids
  stops
}


#' Build CDS lookup vectors from a cds data frame
#'
#' @param cds Data frame with columns: isoform_id, cds_start, cds_stop, strand, coding_status
#' @return Named list with: coding_ids, cds_start_lookup, cds_stop_lookup,
#'   strand_lookup, stop_3prime_lookup
build_cds_lookups <- function(cds) {
  list(
    coding_ids = cds$isoform_id[cds$coding_status == "coding"],
    cds_start_lookup = setNames(cds$cds_start, cds$isoform_id),
    cds_stop_lookup = setNames(cds$cds_stop, cds$isoform_id),
    strand_lookup = setNames(cds$strand, cds$isoform_id),
    stop_3prime_lookup = setNames(
      ifelse(cds$strand == "-", cds$cds_start, cds$cds_stop),
      cds$isoform_id)
  )
}


# ==============================================================================
# Profile extraction helpers
# ==============================================================================

#' Extract all detailed_events from profiles into a single data frame
#'
#' @param profiles Data frame with a detailed_events list-column
#' @param cols Character vector of columns to keep (NULL = all)
#' @return Data frame with all events concatenated
extract_all_events <- function(profiles, cols = NULL) {
  events_list <- lapply(seq_len(nrow(profiles)), function(i) {
    de <- profiles$detailed_events[[i]]
    if (is.null(de) || nrow(de) == 0) return(NULL)
    if (!is.null(cols)) de <- de[, intersect(cols, names(de)), drop = FALSE]
    de
  })
  do.call(rbind, events_list)
}
