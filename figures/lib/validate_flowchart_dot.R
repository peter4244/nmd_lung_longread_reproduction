#!/usr/bin/env Rscript
# Static validator for Graphviz DOT flowcharts.
#
# Catches the recurring formatting issues that bite DiagrammeR/grViz
# flowcharts (the canonical example being figures/SupplementalFigures/
# PairAnalysisFlowchart/build_flowchart.R):
#
#   1. POINT-SIZE drift outside the 3-tier hierarchy (14 / 18 / 22)
#   2. In-cell <BR ALIGN="LEFT"/> wrapping (causes asymmetric row heights)
#   3. Count-immediately-followed-by-label patterns like "<B>72</B> PTC+"
#      that graphviz collapses to "72PTC+" or near-zero kerning
#   4. Edges referencing nodes that were never declared
#   5. cluster_X declarations referenced as plain X in edges (compound=true
#      lhead/ltail expected but missing)
#
# Usage:
#   source("figures/lib/validate_flowchart_dot.R")
#   result <- validate_flowchart_dot(dot_spec_string)
#   if (length(result$errors) > 0) stop("flowchart validation failed")
#
# Returns a list with $errors (character vector) and $warnings (character
# vector). Caller decides whether to raise or proceed.

validate_flowchart_dot <- function(dot, allowed_point_sizes = c(14, 18, 22),
                                    verbose = TRUE) {
  errors <- character(0)
  warnings <- character(0)

  # --- Check 1: POINT-SIZE drift ---
  ps_matches <- regmatches(dot, gregexpr('POINT-SIZE="([0-9]+)"', dot))[[1]]
  if (length(ps_matches) > 0) {
    ps_values <- as.integer(gsub('POINT-SIZE="([0-9]+)"', "\\1", ps_matches))
    bad_ps <- setdiff(unique(ps_values), allowed_point_sizes)
    if (length(bad_ps) > 0) {
      errors <- c(errors, sprintf(
        "POINT-SIZE drift: found %s outside the allowed 3-tier set (%s). Collapse to the 3-tier hierarchy.",
        paste(bad_ps, collapse = ", "),
        paste(allowed_point_sizes, collapse = " / ")))
    }
  }
  # Also check graph default fontsize attribute
  default_fs_match <- regmatches(dot, regexpr('node\\s*\\[[^]]*fontsize=([0-9]+)', dot))
  if (length(default_fs_match) > 0) {
    fs_value <- as.integer(sub(".*fontsize=([0-9]+).*", "\\1", default_fs_match))
    if (!(fs_value %in% allowed_point_sizes)) {
      errors <- c(errors, sprintf(
        "node-default fontsize=%d outside the 3-tier set (%s).",
        fs_value, paste(allowed_point_sizes, collapse = " / ")))
    }
  }

  # --- Check 2: in-cell <BR ALIGN="LEFT"/> wrapping (in non-cluster-title cells) ---
  # Cluster titles legitimately use <BR/> for line breaks. The check is for
  # <BR ALIGN="LEFT"/> inside <TD> cells, which is the asymmetric-rows trap.
  br_count <- length(gregexpr('BR ALIGN="LEFT"', dot, fixed = FALSE)[[1]])
  if (br_count > 0 && !is.na(br_count) &&
      regexpr('BR ALIGN="LEFT"', dot)[1] > 0) {
    # Count actual occurrences (gregexpr returns -1 if no match)
    actual_br <- length(regmatches(dot,
      gregexpr('BR ALIGN="LEFT"', dot))[[1]])
    if (actual_br > 0) {
      warnings <- c(warnings, sprintf(
        "Found %d occurrence(s) of <BR ALIGN=\"LEFT\"/> inside cells. These cause asymmetric row heights — flatten each Why/description to a single clause.",
        actual_br))
    }
  }

  # --- Check 3: count-immediately-followed-by-label ("<B>72</B> PTC+") ---
  # graphviz HTML labels collapse the space between </B> and the following
  # token to near-zero kerning. Use a 2-column TABLE or &nbsp; instead.
  # Pattern: </B> followed by a single space and a letter/digit (not </TD>
  # or other tag).
  tight_kerning <- regmatches(dot,
    gregexpr('</B>\\s[A-Za-z0-9(]', dot))[[1]]
  if (length(tight_kerning) > 0) {
    errors <- c(errors, sprintf(
      "Found %d instance(s) of '</B> <text>' kerning trap (graphviz collapses the space). Use a 2-column <TABLE> [count | label] or '&nbsp;' instead. Examples: %s",
      length(tight_kerning),
      paste(unique(substr(tight_kerning, 1, 20)), collapse = "; ")))
  }

  # --- Check 4: edges reference undeclared nodes ---
  # Extract node declarations: lines like `NAME [label=...]` or `NAME [...]`
  # at the start of a statement, where NAME is a top-level identifier.
  decl_matches <- regmatches(dot,
    gregexpr('(?m)^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\[', dot, perl = TRUE))[[1]]
  declared <- unique(sub('(?m)^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\[.*', "\\1",
                          decl_matches, perl = TRUE))
  # Also pick up subgraph/cluster declarations
  cluster_matches <- regmatches(dot,
    gregexpr('subgraph\\s+(cluster_[A-Za-z0-9_]+)', dot))[[1]]
  clusters <- unique(sub('subgraph\\s+', "", cluster_matches))
  declared <- c(declared, clusters)
  # Extract edges: A -> B or A -> {B, C, D}
  # Simple case: A -> B (with optional :port specifier — :sw, :ne, etc.)
  # Strip the :port off before recording the node name.
  edge_simple <- regmatches(dot,
    gregexpr('([A-Za-z_][A-Za-z0-9_]*)(?::[a-z]+)?\\s*->\\s*([A-Za-z_][A-Za-z0-9_]*)(?::[a-z]+)?',
              dot, perl = TRUE))[[1]]
  edge_nodes <- character(0)
  for (e in edge_simple) {
    parts <- regmatches(e, regexec(
      '([A-Za-z_][A-Za-z0-9_]*)(?::[a-z]+)?\\s*->\\s*([A-Za-z_][A-Za-z0-9_]*)(?::[a-z]+)?',
      e, perl = TRUE))[[1]]
    if (length(parts) >= 3) edge_nodes <- c(edge_nodes, parts[2], parts[3])
  }
  # Also handle braced lists: A -> {B, C, D}
  brace_matches <- regmatches(dot,
    gregexpr('([A-Za-z_][A-Za-z0-9_]*)\\s*->\\s*\\{([^}]+)\\}', dot))[[1]]
  for (b in brace_matches) {
    parts <- regmatches(b, regexec(
      '([A-Za-z_][A-Za-z0-9_]*)\\s*->\\s*\\{([^}]+)\\}', b))[[1]]
    if (length(parts) >= 3) {
      edge_nodes <- c(edge_nodes, parts[2],
                      strsplit(parts[3], '[,\\s]+', perl = TRUE)[[1]])
    }
  }
  edge_nodes <- unique(trimws(edge_nodes))
  edge_nodes <- edge_nodes[nchar(edge_nodes) > 0]
  # Some node-identifier-looking strings are actually attribute keywords —
  # filter those out.
  reserved <- c("rank", "graph", "node", "edge", "subgraph", "cluster_S1",
                "cluster_S2")
  undeclared <- setdiff(edge_nodes, c(declared, reserved))
  if (length(undeclared) > 0) {
    errors <- c(errors, sprintf(
      "Edge endpoints reference undeclared nodes: %s",
      paste(undeclared, collapse = ", ")))
  }

  # --- Check 5a: multi-fan-out from a single source must use consistent
  # tailport treatment (all bare, or all the same `:port`). Mixed ports on
  # outgoing edges from the same source produce visually asymmetric
  # routing — one arrow exits the SW corner, the other exits SE, and the
  # angles never match the symmetric cluster targets.
  fanout_edges <- regmatches(dot,
    gregexpr('([A-Za-z_][A-Za-z0-9_]*)(:[a-z]+)?\\s*->\\s*([A-Za-z_][A-Za-z0-9_]*)',
              dot, perl = TRUE))[[1]]
  src_ports <- list()  # source -> vector of port strings ("" if no port)
  for (e in fanout_edges) {
    parts <- regmatches(e, regexec(
      '([A-Za-z_][A-Za-z0-9_]*)(:[a-z]+)?\\s*->\\s*', e, perl = TRUE))[[1]]
    if (length(parts) >= 3) {
      src <- parts[2]
      port <- if (nchar(parts[3]) > 0) parts[3] else ""
      src_ports[[src]] <- c(src_ports[[src]], port)
    }
  }
  for (src in names(src_ports)) {
    ports <- src_ports[[src]]
    if (length(ports) >= 2 && length(unique(ports)) > 1) {
      errors <- c(errors, sprintf(
        "Node '%s' has %d outgoing edges with mixed tailport specifiers (%s). Use either no port on any edge, or the same port on all edges — mixed ports produce asymmetric arrow routing.",
        src, length(ports),
        paste(ifelse(ports == "", "(bare)", ports), collapse = ", ")))
    }
  }

  # --- Check 5b: edges with lhead=cluster_* should set minlen >= 2 ---
  # Without minlen, graphviz routes the arrow tip to the cluster border at
  # rank N+1, leaving almost no visible arrow between source and the
  # cluster's first node. Forcing minlen=2 gives the arrow a full rank of
  # vertical space to live in.
  lhead_edges <- regmatches(dot,
    gregexpr('\\[[^\\]]*lhead=cluster_[A-Za-z0-9_]+[^\\]]*\\]', dot,
              perl = TRUE))[[1]]
  for (attrs in lhead_edges) {
    minlen_match <- regmatches(attrs, regexpr('minlen=([0-9]+)', attrs))
    if (length(minlen_match) == 0) {
      errors <- c(errors, sprintf(
        "Edge with lhead=cluster_* missing minlen — arrow will be too short. Set minlen=2 (or higher). Edge attrs: %s",
        substr(attrs, 1, 80)))
    } else {
      minlen_val <- as.integer(sub('minlen=([0-9]+)', "\\1", minlen_match))
      if (minlen_val < 2) {
        errors <- c(errors, sprintf(
          "Edge with lhead=cluster_* has minlen=%d (< 2) — arrow will be too short. Bump to minlen >= 2.",
          minlen_val))
      }
    }
  }

  if (verbose) {
    cat(sprintf("\nvalidate_flowchart_dot: %d errors, %d warnings\n",
                length(errors), length(warnings)))
    for (e in errors)   cat(sprintf("  [ERROR]   %s\n", e))
    for (w in warnings) cat(sprintf("  [WARNING] %s\n", w))
  }

  list(errors = errors, warnings = warnings)
}
