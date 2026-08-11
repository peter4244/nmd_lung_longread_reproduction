# =============================================================================
# load_config.R — resolve input roots without hardcoding paths.
#
#   source("R/load_config.R")
#   P <- nmd_paths()
#   readRDS(file.path(P$LEGACY, "nmd_fig_data", "dge_isoform_longread_2026.3.3.rds"))
#
# Resolution per key: environment variable > config/paths.yml > repo-relative.
# =============================================================================

nmd_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  while (!file.exists(file.path(d, "config", "paths.yml"))) {
    p <- dirname(d)
    if (identical(p, d)) stop("cannot locate repository root (no config/paths.yml above ", start, ")")
    d <- p
  }
  d
}

nmd_paths <- function(root = nmd_repo_root()) {
  f <- file.path(root, "config", "paths.yml")
  ln <- grep("^[A-Z_]+:", readLines(f, warn = FALSE), value = TRUE)
  keys <- sub(":.*$", "", ln)
  vals <- trimws(sub("#.*$", "", sub("^[A-Z_]+: *", "", ln)))
  out <- stats::setNames(as.list(vals), keys)
  raw <- out                       # the declared values, before any override
  from_env <- character(0)
  for (k in names(out)) {
    ev <- Sys.getenv(paste0("NMD_", k), unset = NA)
    used_env <- !is.na(ev) && nzchar(ev)
    if (used_env) from_env <- c(from_env, k)
    v  <- if (used_env) ev else out[[k]]
    out[[k]] <- if (grepl("^(/|~)", v)) path.expand(v) else normalizePath(file.path(root, v), mustWork = FALSE)
  }

  # WHY THE KEYS ARE **NOT** COUPLED, having tried it and broken the pipeline.
  #
  # The plan review's A2 said a nested output root should follow the root it is declared under,
  # because `NMD_OUT=<scratch>` moves OUT and leaves ISOPAIR_OUT pinned inside the repo -- a
  # "projection" in which 05t and 05u read the warm tree and exit 0. That hazard is real.
  #
  # The coupling is the WRONG FIX and I shipped it before checking. 05s_orfik_scan.R:48-51 says
  # outright that ".P$ISOPAIR_OUT is the declared intermediates root, NOT .P$OUT" and that the two
  # are "separately overridable (NMD_ISOPAIR_OUT / NMD_OUT) for exactly that reason":
  # 05s_b_orfik_scan_extend.R runs 05s as a child with NMD_OUT redirected to a scratch tree so the
  # subset scan lands there, while the child must still read the REAL nmd_classification.rds,
  # cds.rds and structures.rds. Coupling made the child look for those in an empty directory and
  # the extension died with "cannot open the connection" -- measured 2026-07-29.
  #
  # So the separation stays, and the hazard is caught by nmd_assert_projection() instead: a LOUD
  # check at the top of a run that claims to be a projection, rather than a silent change to
  # semantics a script deliberately depends on. A guard is the right instrument for "you did not
  # mean this"; redefining resolution is not.
  out$ROOT <- root
  attr(out, "from_env") <- from_env
  out
}

# =============================================================================
# nmd_assert_projection(.P) — refuse to run a "clean-tree" run that is not one.
#
# The plan review's highest-risk finding: Phase 1's projection fails GREEN. This is the
# check that makes it fail LOUD instead. Call it at the top of any run that claims to be
# executing into a throwaway tree.
#
# It asserts the property the projection is FOR: no declared output root still resolves
# inside the repository once OUT has been moved out of it. A run that satisfies this cannot
# be silently reading the warm tree, which is the only thing the projection buys.
# =============================================================================
nmd_assert_projection <- function(.P, quiet = FALSE) {
  root <- normalizePath(.P$ROOT, mustWork = FALSE)
  inside <- function(p) {
    p <- normalizePath(p, mustWork = FALSE)
    identical(p, root) || startsWith(p, paste0(root, "/"))
  }
  if (inside(.P$OUT))
    stop("nmd_assert_projection: OUT still resolves inside the repository (", .P$OUT,
         "). This is not a projection -- set NMD_OUT to a directory outside ", root, ".")
  # DATA roots only. FIGURES is deliberately excluded: figures/ is tracked source that SHIPS
  # (D1, amended to admit figures to the citable repo), so a projection legitimately writes
  # renders back into the repo. Requiring it outside would make the guard unsatisfiable and it
  # would be disabled -- which is how guards come to pass while the defect they name is live.
  OUTPUT_ROOTS <- intersect(c("OUT", "CACHE", "ISOPAIR_OUT"), names(.P))
  bad <- Filter(function(k) inside(.P[[k]]), OUTPUT_ROOTS)
  if (length(bad))
    stop("nmd_assert_projection: OUT was redirected but these DATA roots are still inside the ",
         "repository: ", paste(sprintf("%s=%s", bad, unlist(.P[bad])), collapse = ", "),
         ". The run would read the warm tree and exit 0.")
  if (!quiet) cat("[projection] OK -- all output roots resolve outside", root, "\n")
  invisible(TRUE)
}

# =============================================================================
# nmd_assert_writable(dir, .P) — refuse to write into an input root.
#
# WHY THIS IS CODE AND NOT A LINT RULE. tools/lint_paths.sh tries to prove the same
# thing statically, by following assignments from .P$OUT through one level of
# indirection. That cannot work once a script takes `--data-dir`: the destination is
# not knowable until runtime, and the honest static answer is "flag it", which is why
# that check had been red for weeks and was therefore being ignored.
#
# The rule is the same either way -- a pipeline must never write into a tree it reads
# as input -- so put it where it can actually be decided.
#
# Compares on PATH-COMPONENT boundaries, deliberately: "/data" is a string prefix of
# "/data_v2", and a prefix test would silently pass a write into a sibling tree.
# =============================================================================
nmd_assert_writable <- function(dir, .P, what = "output directory") {
  INPUT_ROOTS <- c("DEPOSIT", "SQANTI", "ISOCALL", "ANNOT", "REFERENCE", "LEGACY",
                   "FIGDATA", "MASHR_SR", "MASHR_LR", "ISOPAIR", "PHENO", "MODEL")
  # RESOLVE THROUGH THE NEAREST EXISTING ANCESTOR.
  #
  # normalizePath(mustWork = FALSE) does NOT resolve symlinks through a path that does not
  # exist yet -- and the destination usually does not exist yet, because the caller is about
  # to dir.create() it. This repo reaches an input root through a symlink (data_deposit ->
  # nmd_deposit_2026), so the guard compared an UNRESOLVED alias against a RESOLVED root,
  # found no shared prefix, and allowed the write. Found by review 2026-07-26; the earlier
  # "verified both directions" test only exercised an existing path via the canonical route,
  # which is the one case that could not fail.
  d <- normalizePath(dir, mustWork = FALSE)
  if (!dir.exists(dir)) {
    anc <- dir; tail <- character(0)
    while (!dir.exists(anc) && nzchar(anc) && anc != dirname(anc)) {
      tail <- c(basename(anc), tail); anc <- dirname(anc)
    }
    if (dir.exists(anc)) d <- do.call(file.path, c(list(normalizePath(anc)), as.list(tail)))
  }
  for (k in INPUT_ROOTS) {
    r <- .P[[k]]
    if (is.null(r) || !nzchar(r)) next
    r <- normalizePath(r, mustWork = FALSE)
    # equal, or d is strictly below r on a component boundary
    if (identical(d, r) || startsWith(d, paste0(r, .Platform$file.sep))) {
      stop(sprintf(paste0(
        "%s resolves INSIDE the input root .P$%s\n",
        "    %s\n  is under\n    %s\n",
        "  Writing there overwrites the reference copies every verification compares against.\n",
        "  Point it at .P$OUT (or a directory beneath it) instead."), what, k, d, r),
        call. = FALSE)
    }
  }
  invisible(d)
}

# =============================================================================
# nmd_data_dir(.P) — where the isopair intermediates live, for anything that READS them.
#
# WHY THIS IS A FUNCTION AND NOT TWENTY COPIED LINES. Ten figure exporters each resolved
# this themselves, and every one of them defaulted to .P$ISOPAIR -- the LEGACY tree. That is
# not a style problem: the legacy data_mashr carries the 2026-03-11 structures.rds, whose
# isoform set differs from the rebuilt one, so those scripts exported panels from a
# different universe than the reports use. The paper's main text and its supplement could
# disagree and nothing would say so. W10 fixed one of them by hand; the remaining nine would
# have drifted again the moment anyone copied a tenth.
#
# Default .P$OUT/data_mashr, overridable with --data-dir for stage-by-stage verification
# against retained references. Fails loudly when absent, because a missing intermediate must
# stop the run rather than silently fall back to a copy from another vintage -- that silent
# fallback is the defect this whole exercise exists to document.
# =============================================================================

# =============================================================================
# nmd_ensure_out(.P) -- create the DECLARED OUTPUT roots.
#
# config/paths.yml says of CACHE and OUT: "outputs and are created on demand". NOTHING
# IMPLEMENTED THAT. The authoring tree has had a tmp/out since the first run, so every script
# writing into it worked here and nowhere else -- and the clean-room rerun found it the
# expensive way on 2026-07-28: Isoform_Level_Quantification.Rmd:301 does
# `saveRDS(dge_nofilt, file.path(.P$OUT, ...))` and the only dir.create that would have made
# .P$OUT sits at :593, 290 lines later. It failed after 376 seconds of real work with
# `gzfile(): cannot open the connection`, and took the whole §1-§2 chain down behind it: two
# steps failed writing, four more failed reading what those two never wrote.
#
# ONE PLACE RATHER THAN A dir.create IN EVERY WRITER, because the next script to write into
# .P$OUT would reintroduce it. The roots here are OUTPUTS by declaration -- nmd_assert_writable
# already refuses any write that lands inside an INPUT root -- so creating them is exactly what
# the config's own contract promises.
#
# Not called from nmd_paths(): resolving a path must not touch the filesystem, or read-only
# tooling (tools/project.py's gate, tools/dryrun.py) starts creating directories in the tree it
# is auditing.
# =============================================================================
nmd_ensure_out <- function(.P, extra = character(0)) {
  for (d in c(.P$OUT, .P$CACHE, extra)) {
    if (!is.null(d) && nzchar(d) && !dir.exists(d))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(.P)
}

nmd_data_dir <- function(.P, what = NULL) {
  # .P$ISOPAIR_OUT when declared (config/paths.yml), else .P$OUT/data_mashr. The default
  # used to be data_mashr unconditionally, which on 2026-07-26 silently resolved to the
  # UNFLOORED tree and produced pop_BC = 2,942 against an expected ~1,400-1,550.
  d <- if (!is.null(what)) file.path(.P$OUT, what)
       else if (!is.null(.P$ISOPAIR_OUT) && nzchar(.P$ISOPAIR_OUT)) .P$ISOPAIR_OUT
       else file.path(.P$OUT, "data_mashr")
  a <- commandArgs(trailingOnly = TRUE)
  if ("--data-dir" %in% a) {
    i <- which(a == "--data-dir")
    if (i < length(a)) d <- a[i + 1]
  }
  d <- normalizePath(d, mustWork = FALSE)
  if (!dir.exists(d)) {
    stop(sprintf(paste0(
      "isopair data directory does not exist:\n    %s\n",
      "  Run 01_prepare_data_mashr.R -> 01b_build_isoform_infrastructure.R -> ",
      "02_build_profiles_mashr.R first,\n  or pass --data-dir <path> to read a retained tree."),
      d), call. = FALSE)
  }
  # ANNOUNCE THE POPULATION. A wrong tree is invisible until a panel number is already wrong;
  # printing the C4 pair count makes it visible at resolution time, for free. On 2026-07-26 an
  # unfloored tree resolved silently and gave pop_BC 2,942 against an expected ~1,400-1,550.
  pf <- file.path(d, "profiles_c4_allsamples.rds")
  if (file.exists(pf)) {
    n <- tryCatch(nrow(readRDS(pf)), error = function(e) NA_integer_)
    cat(sprintf("isopair data dir: %s  (C4 pairs: %s)\n", d,
                if (is.na(n)) "unreadable" else format(n, big.mark = ",")))
  }
  d
}

# =============================================================================
# Companion SSOT helpers, sourced here so that `source("R/load_config.R")` -- which every
# analysis script already does -- is sufficient.
#
# R/mashr_columns.R defines nmd_mash_cols(). It was created to prevent a SILENT ZERO: the
# mashr objects carry two cell-type spellings, FB and MV are identical in both, so a
# mismatch resolves two of four cell types and silently halves a result (2026-07-25: a
# reported range collapsed from 20.7-54.5% to 0-38.9% with no warning).
#
# Nothing sourced it. Three production scripts called nmd_mash_cols() anyway --
# predictor_comparison/01, isopair_wrapper/build_paper_tables.R, and
# 05_final_report_mashr.Rmd -- and every one of them died with "could not find function".
# Found 2026-07-26 by running 01. A guard nothing loads is not a guard.
# =============================================================================
local({
  f <- file.path(nmd_repo_root(), "R", "mashr_columns.R")
  # FAIL LOUDLY IF IT IS ABSENT (2026-08-03). This was `if (file.exists(f)) source(f)`, which
  # made absence SILENT -- and absence is not benign here. This file is the SSOT for resolving
  # mashr cell-type columns across two naming conventions; without it FB and MV still resolve
  # while AT2 and LAE do not, so the run reports two of four cell types and does not say so. The
  # header of mashr_columns.R records the cost: a published range moved from 20.7-54.5% to
  # 0-38.9% with no warning. A clean-room projection built from a keep-list that omitted this
  # file would have reproduced that silently (C93), which is exactly what D14 forecloses -- a
  # missing file must surface as a crash, not as a quietly wrong number.
  if (!file.exists(f)) {
    stop("R/mashr_columns.R not found at ", f, ".\n",
         "It resolves mashr cell-type column names and its absence does NOT fail loudly on its ",
         "own: two of four cell types would silently drop out. Restore it before running.")
  }
  # INLINE, NOT source(f). extract detects source() with a CONSTRUCTED PATH argument but not
  # source(VARIABLE, local = ...) -- so written this way the dependency appears in the graph
  # and mashr_columns.R is projected into a clean room instead of silently missing (W277).
  source(file.path(nmd_repo_root(), "R", "mashr_columns.R"), local = FALSE)
})
