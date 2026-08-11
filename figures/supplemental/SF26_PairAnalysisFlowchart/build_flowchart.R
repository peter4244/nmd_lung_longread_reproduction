#!/usr/bin/env Rscript
# Pair-analysis cohort flowchart — standalone supplement render.
#
# Plain-language flow chart of how the pair analysis arrives at the gene-matched
# NMD-susceptible / Control pairs (upstream chain) and how those pairs
# then narrow into the two analysis subsets used in the manuscript
# (downstream cascade), with kept / dropped counts and the reason for each
# dropout at every filter step.
#
# Outputs (in this folder):
#   figure_s_pair_analysis_flowchart.dot   — raw Graphviz DOT spec
#   figure_s_pair_analysis_flowchart.html  — standalone htmlwidget; opens
#                                            in any browser
#
# For PNG / PDF: open the HTML in Chrome and File → Print → "Save as PDF",
# or paste the .dot file into https://dreampuf.github.io/GraphvizOnline/.
#
# EVERY number is derived by cascade_counts.R into data/cascade_counts.tsv; this script
# only formats them. Run cascade_counts.R first, or this stops.
#
# The old header pointed at 05_final_report_gencode_scope_2026-07-11.Rmd (a path that no
# longer exists) and named two verifiers as the drift guard. Those verifiers hold the
# SUPERSEDED values as their expected manifest, so they did not merely fail to catch the
# 2026-08-08 drift -- they certified it. Do not treat them as the guard here. The guard is
# that nothing is typed: cascade_counts.R asserts all 13 cascade sums close, and N() stops
# on a key it does not have.

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(DiagrammeR)
  library(htmlwidgets)
})

if (Sys.getenv("RSTUDIO_PANDOC") == "" && file.exists(
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc")) {
  Sys.setenv(RSTUDIO_PANDOC =
             "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
}

# Repointed 2026-07-26 (W31 step 2). Was NMD_ROOT/results/.../data_mashr under .P$LEGACY,
# so the flowchart counted the PUBLISHED pair set regardless of what had been rebuilt.
DM <- nmd_data_dir(.P)
HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- args[grep("--file=", args)]
  if (length(m) > 0) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  normalizePath(getwd())
})()

# ─────────────────────────────────────────────────────────────────────
# EVERY count comes from cascade_counts.R. None is typed here.
#
# 2026-08-09: until today the upstream nodes read nrow(profiles_c2) live while every
# subset box carried a hardcoded string, so the 2026-08-08 deposit rebuild updated half
# the figure and left the rest. The rendered .dot printed 1,548 TWICE and 1,578 ONCE, with
# the 1,578 arrow pointing into the box labelled 1,548, and one box gave the
# reference-share floor as 5,001 three lines below reporting 5,049 candidates.
#
# Retyping the literals would rebuild the same trap on the next substrate change, so they
# are gone. `N()` stops on an unknown key rather than returning NA, which would render as
# an empty box -- a missing number must be impossible to ship, not merely ugly.
# ─────────────────────────────────────────────────────────────────────
CASCADE_TSV <- file.path(HERE, "data", "cascade_counts.tsv")
if (!file.exists(CASCADE_TSV))
  stop("missing ", CASCADE_TSV, "\n  Run: Rscript cascade_counts.R")
cascade <- read.delim(CASCADE_TSV, stringsAsFactors = FALSE)
N <- function(key) {
  i <- match(key, cascade$key)
  if (is.na(i)) stop(sprintf("cascade_counts.tsv has no key '%s' -- rerun cascade_counts.R", key))
  as.integer(cascade$value[i])
}
fmt  <- function(x) format(x, big.mark = ",")
fmtN <- function(key) fmt(N(key))

N_RAW_ISO    <- N("n_raw_iso")
N_FILTER_ISO <- N("n_filter_iso")
N_4CT_SMP    <- N("n_4ct_samples")
N_NMD_AS     <- N("n_nmd_as")
N_NONNMD_AS  <- N("n_nonnmd_as")
N_C2         <- N("n_c2")
N_C4         <- N("n_c4")

# ─────────────────────────────────────────────────────────────────────
# DOT spec — three logical zones:
#   1. Upstream pipeline (single column, top)
#   2. Matched-pair hub (pop_BC, both arms)
#   3. Two subset clusters (left = strict / right = broader), each with
#      cluster-level header as its title and filter nodes inside
#   + Hidden-PTC sub-cascade off Subset 2
#
# Box formatting: each filter node uses a small HTML table inside the
# label so the kept / dropped / why rows align cleanly.
# ─────────────────────────────────────────────────────────────────────
filter_node <- function(what, kept, dropped, why) {
  # 2-column TABLE — left column = label, right column = value.
  # Eliminates the "<B>72</B> PTC+" kerning trap (graphviz collapses the
  # space after </B>); proper grid layout aligns counts visually.
  sprintf(
    '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Filter:</B>&nbsp;%s</TD></TR>
       <TR>
         <TD ALIGN="LEFT"><FONT COLOR="#1a5a2a"><B>kept:</B></FONT></TD>
         <TD ALIGN="LEFT"><FONT COLOR="#1a5a2a"><B>%s</B></FONT></TD>
       </TR>
       <TR>
         <TD ALIGN="LEFT"><FONT COLOR="#b22222">dropped:</FONT></TD>
         <TD ALIGN="LEFT"><FONT COLOR="#b22222">%s</FONT></TD>
       </TR>
       <TR>
         <TD ALIGN="LEFT" COLSPAN="2"><FONT POINT-SIZE="18" COLOR="#444">Why: %s</FONT></TD>
       </TR>
     </TABLE>>', what, kept, dropped, why)
}

upstream_node <- function(stage, detail, what_dropped = NULL, center = FALSE) {
  al <- if (center) "CENTER" else "LEFT"
  rows <- sprintf(
    '<TR><TD ALIGN="%s"><B>%s</B></TD></TR><TR><TD ALIGN="%s">%s</TD></TR>',
    al, stage, al, detail)
  if (!is.null(what_dropped)) {
    rows <- paste0(rows, sprintf(
      '<TR><TD ALIGN="%s"><FONT POINT-SIZE="18" COLOR="#b22222">%s</FONT></TD></TR>',
      al, what_dropped))
  }
  sprintf('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">%s</TABLE>>', rows)
}

# HUB_CONS lists Figure 3 Panels B and C as of 2026-08-09. They were missing, and both run on
# pop_BC with no CDS restriction -- Panel B is 1,578 per arm, and Panel C's denominators are 1,578
# (SE 809/1,578 = 51.3% against 342/1,578 = 21.7%). The omission mirrors the published Figure 3
# legend, which wrongly restricts Panel C to the GENCODE subset (claim 4.4.2).
#
# KEEP EXPLANATIONS OUT OF THE TEMPLATE BELOW. R's sprintf caps `fmt` at 8192 characters and this
# string is close to it; four lines of DOT comment tipped it over and the failure reads as
# "'fmt' length exceeds maximal format length 8192", printed underneath a full dump of the
# template, which buries it. Two other traps in the same string, both hit today: a literal
# apostrophe ends the single-quoted R string, and a literal % is read as a format specifier.

dot_spec <- sprintf('
digraph cohort_flow {
  graph [rankdir=TB, fontname="Helvetica", nodesep=0.20, ranksep=0.35,
         splines=ortho, compound=true]
  node  [shape=box, style="rounded,filled", fontname="Helvetica",
         fontsize=18, margin="0.15,0.12"]
  edge  [color="#666", penwidth=1.4, fontname="Helvetica",
         fontsize=18, fontcolor="#444"]

  // ─── Starting cohort + upstream classification + pair construction
  U1 [label=%s, fillcolor="#e8eef4", color="#2f4858"]
  U2 [label=%s, fillcolor="#fff7bc", color="#d95f0e"]
  U3 [label=%s, fillcolor="#fff7bc", color="#d95f0e"]

  U1 -> U2 -> U3

  // ─── Hub: matched isoform pairs ──────────────────────────────────
  HUB [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                <TR><TD ALIGN="CENTER"><FONT COLOR="white" POINT-SIZE="22"><B>Matched isoform pairs</B></FONT></TD></TR>
                <TR><TD ALIGN="CENTER"><FONT COLOR="white"><B>%s NMD susceptible  ·  %s Control</B></FONT></TD></TR>
                <TR><TD ALIGN="CENTER"><FONT COLOR="#cfd8e2" POINT-SIZE="18">one triplet per gene</FONT></TD></TR>
              </TABLE>>,
       fillcolor="#1f3a5f", color="#0a1a30", penwidth=2.0]

  U3 -> HUB

  // ─── Subset 1 cluster (left) ─────────────────────────────────────
  subgraph cluster_S1 {
    label = <<FONT POINT-SIZE="22"><B>GENCODE-restricted subset (n = %s)</B></FONT><BR/><BR/>all three transcripts in the triplet are annotated coding transcripts>
    labelloc = "t"
    style = "rounded,filled"
    fillcolor = "#fef0e6"
    color = "#7d4030"
    fontcolor = "#5a2d1e"
    margin = 12

    S1_F1 [label=%s, fillcolor="white", color="#a2604a"]
    S1_F2 [label=%s, fillcolor="white", color="#a2604a"]

    S1_GROUPS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Classification:</B>&nbsp;premature stop in the NMD-susceptible transcript?</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">PTC+ (%s%%)</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">PTC−</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">Control</TD></TR>
                     </TABLE>>,
               fillcolor="#fff7bc", color="#d95f0e"]

    S1_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 3 (Panels D / E / F)</TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 4 (Panels A / B)</TD></TR>
                     <TR><TD ALIGN="LEFT">SF33</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

    S1_F1 -> S1_F2 -> S1_GROUPS -> S1_CONS [style=dashed, color="#7d4030"]
  }

  // ─── Subset 2 cluster (right) ────────────────────────────────────
  subgraph cluster_S2 {
    label = <<FONT POINT-SIZE="22"><B>Reference AUG-traceable subset (n = %s)</B></FONT><BR/><BR/>reference is annotated; NMD / Control can be novel if the reference start codon maps to the comparator>
    labelloc = "t"
    style = "rounded,filled"
    fillcolor = "#e6f1f8"
    color = "#1f5570"
    fontcolor = "#0e2c3e"
    margin = 12

    S2_F1 [label=%s, fillcolor="white", color="#406080"]
    S2_F2 [label=%s, fillcolor="white", color="#406080"]
    S2_F3 [label=%s, fillcolor="white", color="#406080"]

    S2_GROUPS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Classification:</B>&nbsp;reference start codon projected to a premature stop?</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">PTC+ (%s%%)</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">PTC−</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>%s</B></TD><TD ALIGN="LEFT">Control</TD></TR>
                     </TABLE>>,
               fillcolor="#fff7bc", color="#d95f0e"]

    S2_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 4 (Panels C / D)</TD></TR>
                     <TR><TD ALIGN="LEFT">SF31 / SF32</TD></TR>
                     <TR><TD ALIGN="LEFT">SF34 / SF36</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

    S2_F1 -> S2_F2 -> S2_F3 -> S2_GROUPS -> S2_CONS [style=dashed, color="#1f5570"]
  }

  // ─── Hidden-PTC sub-cascade ─────────────────────────────────────
  HPTC [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                <TR><TD ALIGN="LEFT"><B>Occult-PTC subset (n = %s)</B></TD></TR>
                <TR><TD ALIGN="LEFT">reference-anchored analysis flags a premature stop,</TD></TR>
                <TR><TD ALIGN="LEFT">but the TransDecoder2 CDS caller does not</TD></TR>
                <TR><TD ALIGN="LEFT"><FONT COLOR="#b22222">%s / %s dropped</FONT> — no disagreement</TD></TR>
              </TABLE>>,
        fillcolor="#a48bbf", color="#4a3556", fontcolor="#1a0a25"]

  HPTC_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">SF35</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

  // ─── Direct-from-hub consumers (whole pair set, no further narrowing) ─
  HUB_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                    <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                    <TR><TD ALIGN="LEFT">Figure 3 (Panels B / C)</TD></TR>
                    <TR><TD ALIGN="LEFT">SF27 / SF28 / SF29</TD></TR>
                    <TR><TD ALIGN="LEFT">SF30</TD></TR>
                  </TABLE>>,
            fillcolor="#eeeeee", color="#666"]

  // Invisible spacer to force horizontal breathing room on the HUB→HUB_CONS
  // arrow; without it, rank=same pins HUB_CONS one nodesep (0.20) from HUB
  // and the arrow visually collapses.
  HUB_PAD [style=invis, width=1.4, height=0.01, label=""]

  // ─── Edges from hub into subset clusters and TD2 sub-cascade ─────
  // HUB → cluster arrows: no tailport specifier on either edge (mixed ports
  // like sw vs se on a single source produce visually asymmetric routes;
  // default auto-routing under splines=ortho + symmetric cluster targets
  // gives a balanced fan-out). minlen=2 forces each edge to span 2 ranks
  // for visible length. Both rules enforced by
  // figures/lib/validate_flowchart_dot.R.
  HUB -> S1_F1 [lhead=cluster_S1, color="#7d4030", penwidth=2.0, minlen=2]
  HUB -> S2_F1 [lhead=cluster_S2, color="#1f5570", penwidth=2.0, minlen=2]
  // HUB_CONS placed to the RIGHT of HUB (same rank), not below — the
  // whole-pair-set descriptive analyses are a sibling of the two subset
  // narrowings, so a lateral placement reads more accurately.
  HUB -> HUB_PAD [style=invis]
  HUB_PAD -> HUB_CONS [style=invis]
  // Solid arrow — HUB_CONS is a peer consumer of the HUB (same class as
  // HUB → S1 / HUB → S2), so it takes the main-flow edge style. Dashed
  // is reserved for edges inside a cluster (…_GROUPS → …_CONS).
  HUB -> HUB_CONS [color="#555", constraint=false]
  { rank = same; HUB; HUB_PAD; HUB_CONS }

  S2_F3 -> HPTC [color="#4a3556", minlen=3]
  HPTC -> HPTC_CONS [style=dashed, color="#4a3556"]
}',
  upstream_node("Starting cohort",
                sprintf("AT (n = 6) · DD (n = 8) · FB (n = 6) · MV (n = 6) · %d samples<BR/>%s non-fusion isoforms (≥5%% of gene CPM in ≥1 DMSO or SMG1i sample)",
                        N_4CT_SMP, fmt(N_FILTER_ISO))),
  upstream_node("Identify NMD susceptible isoforms",
                sprintf("%s NMD susceptible · %s non-NMD",
                        fmt(N_NMD_AS), fmt(N_NONNMD_AS))),
  upstream_node("Construct paired comparisons",
                sprintf("per gene: dominant NMD-susceptible vs reference · two dominant non-NMD (Control + reference)<BR/>→ %s NMD-susceptible pairs · %s non-NMD Control candidate pairs, gene-matched to the hub below",
                        fmt(N_C2), fmt(N_C4)),
                what_dropped = sprintf("reference-share floor: dropped %s genes (%s → %s) whose selected reference was &lt; 25%% of total gene expression",
                                       fmtN("floor_dropped"), fmtN("floor_before"), fmtN("floor_after")),
                center = TRUE),

  # ── Hub ──
  fmtN("pop_bc_c2"), fmtN("pop_bc_c4"),

  # ── Subset 1: header, two filter nodes, group split ──
  # Why = single clause, no in-cell line wrapping.
  fmtN("secA_n"),
  filter_node("all three transcripts are annotated (in GENCODE)",
              sprintf("%s / %s", fmtN("secA_enst_kept"), fmtN("secA_enst_kept")),
              sprintf("%s / %s", fmtN("secA_enst_dropped"), fmtN("secA_enst_dropped")),
              "most NMD-susceptible transcripts are novel"),
  filter_node("all three transcripts have a curated coding sequence",
              sprintf("%s / %s", fmtN("secA_coding_kept"), fmtN("secA_coding_kept")),
              sprintf("%s / %s", fmtN("secA_coding_dropped"), fmtN("secA_coding_dropped")),
              "annotated non-coding transcripts (lncRNA, etc.)"),
  fmtN("secA_ptc_pos"), N("secA_ptc_pct"), fmtN("secA_ptc_neg"), fmtN("secA_n"),

  # ── Subset 2: header, three filter nodes, group split ──
  fmtN("secC_n"),
  filter_node("reference is annotated (NMD and Control can be novel)",
              sprintf("%s / %s", fmtN("secC_enstref_kept"), fmtN("secC_enstref_kept")),
              sprintf("%s / %s", fmtN("secC_enstref_dropped"), fmtN("secC_enstref_dropped")),
              "reference itself is a novel isoform"),
  filter_node("the reference start codon maps onto the comparator transcript",
              sprintf("%s NMD  ·  %s Control",
                      fmtN("secC_traceable_kept_nmd"), fmtN("secC_traceable_kept_ctrl")),
              sprintf("%s NMD  ·  %s Control",
                      fmtN("secC_traceable_dropped_nmd"), fmtN("secC_traceable_dropped_ctrl")),
              "reference exon missing from comparator"),
  filter_node("keep only genes where both NMD and Control side passed",
              sprintf("%s / %s", fmtN("secC_n"), fmtN("secC_n")),
              sprintf("%s NMD  ·  %s Control",
                      fmtN("secC_match_dropped_nmd"), fmtN("secC_match_dropped_ctrl")),
              "reference passed on one arm but not the other"),
  fmtN("secC_ptc_pos"), N("secC_ptc_pct"), fmtN("secC_ptc_neg"), fmtN("secC_n"),

  # ── Occult-PTC sub-cascade ──
  fmtN("occult_n"), fmtN("occult_dropped"), fmtN("secC_n")
)

dot_path  <- file.path(HERE, "figure_s_pair_analysis_flowchart.dot")
html_path <- file.path(HERE, "figure_s_pair_analysis_flowchart.html")

# Static flowchart validator — catches recurring formatting traps before render
# NMD_ROOT was never defined here -- a leftover from the 2026-07-26 repoint (see the comment at
# the DATA line), which removed the variable and missed this one use. The script therefore ran
# to completion and died at the VALIDATOR, so the figure was unbuildable for the least visible
# possible reason: the check that would have caught a bad figure is what crashed.
source(file.path(.nmd_root, "figures/lib/validate_flowchart_dot.R"))
v <- validate_flowchart_dot(dot_spec)
if (length(v$errors) > 0) {
  stop(sprintf("validate_flowchart_dot found %d error(s) — fix before render",
                length(v$errors)))
}

writeLines(dot_spec, dot_path)
cat(sprintf("[dot]  → %s\n", dot_path))

# Explicit container size matching the SVG's natural extent at the chosen
# font sizes — prevents the htmlwidget's small default container from
# scaling the SVG (which makes fonts look microscopic). Width/height in
# pixels; if you bump font sizes or add more nodes, bump these too.
widget <- DiagrammeR::grViz(dot_spec, width = 1800, height = 2200)
htmlwidgets::saveWidget(widget, html_path,
                         selfcontained = TRUE,
                         title = "Pair-analysis cohort flowchart")
cat(sprintf("[html] → %s\n", html_path))

# Reproducible PNG render via graphviz (same engine DiagrammeR uses via viz.js).
# Falls back to the browser workflow below if `dot` is not on PATH.
png_path <- file.path(HERE, "figure_s_pair_analysis_flowchart.png")
if (nzchar(Sys.which("dot"))) {
  status <- system2("dot", c("-Tpng", "-Gdpi=150", shQuote(dot_path),
                             "-o", shQuote(png_path)))
  if (status == 0) cat(sprintf("[png]  → %s\n", png_path))
} else {
  cat("\n'dot' not found — for PNG/PDF open the HTML in a browser and print.\n")
}
