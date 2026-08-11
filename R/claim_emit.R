#!/usr/bin/env Rscript
# Emit a claim's value from the code that computes it. R side of the emission contract.
#
# The Python side is tools/claim_emit.py and writes the same file in the same format; the
# pipeline is half R and half Python and the two must land in one table or the check has to
# reconcile two shapes, which is where a reconciliation bug would live.
#
# WHY. `g2_value` is prose in 66 of 72 ledger rows, so "does the pipeline still produce the
# paper's number?" cannot be answered by comparison. The previous mechanical attempt grepped
# rendered HTML for the number -- a proxy that CANNOT FAIL -- and scored 20 claims reproduced
# including five known to differ; all 25 outcomes were voided (D19). Here the value is written by
# the producer at the moment it is computed, so there is no gap between the number the code made
# and the number recorded.
#
# WHY `population` IS REQUIRED. Unstated populations are the dominant defect class in this
# manuscript, about one per section. The population is known at the call site and nowhere else,
# so it is captured there or it is lost. Omitting it is an error, not a warning.
#
# Usage, one line where the number is computed:
#
#   source(file.path(.nmd_root, "R", "claim_emit.R"))
#   claim_emit("4.1.3", "genes after the GENCODE-scope filter", nrow(kept),
#              n = nrow(kept),
#              population = "c2 pairs whose reference isoform has a GENCODE CDS; FB excluded")
#
# Output: $NMD_CLAIM_VALUES if set, else tmp/claim_values.tsv, appended. Run in a CLEAN tree --
# appending to a stale file silently mixes vintages.
#
# claim_id IS A CALL-SITE LABEL, NOT THE LEDGER'S ID, AND YOU MUST NOT JOIN ON IT. This was
# undocumented until 2026-08-10 and cost a wrong conclusion the first day the table was read: the
# container seat joined by claim_id, found no row for 4.4.5 or 4.4.6, and reported that neither is
# emitted by the chain and that reading rendered HTML was the only evidence available. Both are
# emitted in full. 4.4.5's five quantities (48/132 = 36.36%, 2/132 = 1.52%, 298 Control events) come
# out under 4.4.7 and 4.4.8; 4.4.6's split (26 / 18 / 4, shares 54.17 / 37.5 / 8.33) comes out TWICE,
# under 4.4.10 and 4.6.1. Found by the paper-validation seat, who confirmed the mapping BY PRODUCER
# CHUNK rather than by value -- which is what makes it a mapping and not a coincidence.
#
# Two facts follow, and both run against the intuition the column name creates:
#   * ONE ledger claim's quantities may be emitted under a DIFFERENT id -- whichever call site
#     happens to compute them.
#   * ONE emitted quantity may serve SEVERAL ledger claims: 4.4.10 and 4.6.1 both carry the
#     mechanism shares from chunk sec2b-attribution, 4.4.7 and 4.5.2 both carry the PTC rates from
#     sec2a-ptc.
# So an id join UNDER-REPORTS coverage, and "absent" from such a join means nothing. JUDGE COVERAGE
# ON quantity AND population, and confirm a mapping through producer_file/producer_line.
#
# THE SHAPE IS THE HAZARD, NOT THIS INSTANCE. The same day, a chunk-label probe reported two chunks
# "absent" from a render that had run them -- that html_document emits no chunk labels at all. Two
# label-based searches, two different artifacts, two false absences, both found by the other seat
# rather than the one that ran them. A label is not the thing it names.
#
# THE SIX §5 COLUMNS (ensemble_member, aggregation, n_background, replicate, sd_within,
# sd_between) exist on this side too and are normally left NULL by §1-4. They are here for
# SHAPE, not for use: both emitters append to one file and whichever writes first sets the
# header, and §1-4 is entirely R -- so if this header were the older ten columns, the Python
# side's six fields would land in csv.DictReader's unnamed restkey and the S4 uniqueness key
# would silently lose ensemble_member. See the comment at `header` below.
#
# `producer_line` is BEST-EFFORT and its value depends on how R was invoked. Under
# `Rscript -e 'rmarkdown::render(...)'` -- the only way §1-4 ever runs -- it is `chunk:<label>`.
# Under `Rscript file.R` srcrefs are absent and it is EMPTY; `producer_file` is still correct.

claim_emit <- function(claim_id, quantity, value, published = NULL, n = NULL,
                       population = NULL, run_id = NULL, ensemble_member = NULL,
                       aggregation = NULL, n_background = NULL, replicate = NULL,
                       sd_within = NULL, sd_between = NULL, restriction = NULL) {
  if (is.null(claim_id) || !nzchar(trimws(as.character(claim_id))))
    stop("claim_emit: claim_id is required")
  if (is.null(quantity) || !nzchar(trimws(as.character(quantity))))
    stop("claim_emit: quantity is required -- two numbers in one claim must be distinguishable")
  if (is.null(population) || !nzchar(trimws(as.character(population))))
    stop(sprintf(paste0("claim_emit: population is required for '%s'. State what was counted and ",
                        "what was excluded, in words. Unstated populations are the most common ",
                        "defect in this manuscript, and the call site is the only place the ",
                        "population is known."), claim_id))

  # SAME CLOSED VOCABULARY AS THE PYTHON SIDE, and for the same reason: the column exists to
  # make |mean_k phi_k| and mean_k |phi_k| distinguishable after the fact, so "mean abs",
  # "abs-mean" and "meanOfAbs" drifting into it would leave exactly the ambiguity it removes.
  if (!is.null(aggregation) &&
      !(trimws(as.character(aggregation)) %in% c("abs_of_mean", "mean_of_abs")))
    stop(sprintf(paste0("claim_emit: aggregation must be 'abs_of_mean' or 'mean_of_abs', got ",
                        "'%s'. abs_of_mean = |mean_k phi_k| (the ensemble's attribution); ",
                        "mean_of_abs = mean_k |phi_k| (average magnitude). If you are unsure ",
                        "which the producer computes, read it before emitting -- guessing here ",
                        "silently mislabels the estimand."), as.character(aggregation)))

  root <- tryCatch({
    d <- getwd()
    while (!file.exists(file.path(d, "config", "paths.yml"))) {
      p <- dirname(d); if (identical(p, d)) stop("no root"); d <- p
    }
    d
  }, error = function(e) getwd())

  target <- Sys.getenv("NMD_CLAIM_VALUES", unset = "")
  if (!nzchar(target)) target <- file.path(root, "tmp", "claim_values.tsv")
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)

  # Tabs and newlines would silently reshape the row into a different record.
  clean <- function(x) gsub("[\t\r\n]", " ", trimws(paste0(as.character(x), collapse = " ")))

  # The call site, so a disagreement can be taken straight to the line that produced it.
  # Rscript parses with keep.source=FALSE, so srcrefs are usually ABSENT and the first version of
  # this recorded "unknown" for every row. The FILE is recoverable regardless, from --file= in the
  # raw command line, so take that as the reliable part and treat the line number as best-effort.
  # Every upstream report is rendered by `Rscript -e 'rmarkdown::render(...)'`, which has NO
  # --file=, so the first version recorded "interactive" for all of them -- i.e. for the majority
  # of the pipeline, the producer column of the reader-facing index would have been empty.
  # knitr::current_input() is the reliable answer inside a knit; fall back through the call stack.
  src_file <- local({
    if (requireNamespace("knitr", quietly = TRUE)) {
      ci <- tryCatch(knitr::current_input(), error = function(e) NULL)
      if (!is.null(ci) && nzchar(ci)) return(basename(ci))
    }
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) return(basename(sub("^--file=", "", a[1])))
    for (i in rev(seq_len(sys.nframe()))) {
      sf <- tryCatch(attr(attr(sys.frame(i), "srcref"), "srcfile"), error = function(e) NULL)
      if (!is.null(sf$filename) && nzchar(sf$filename)) return(basename(sf$filename))
    }
    "unknown"
  })
  # THE LOCATOR. A line number if one is recoverable, otherwise the KNITR CHUNK LABEL, prefixed
  # `chunk:` so it can never be misread as a line. Measured 2026-07-29: under the pipeline's real
  # invocation -- `Rscript -e 'rmarkdown::render(...)'`, which is how rerun_roles.tsv renders
  # every §1-4 report -- srcrefs are absent and this column came back EMPTY, so all 57 §1-4
  # emissions would have printed as `Report.Rmd:` with nothing to take to the code. That is the
  # "a disagreement you cannot take to a line" defect C50 refused to accept when it was proposed
  # as a rename, arriving instead by omission.
  #
  # The chunk label is the BETTER address here, not a consolation: docs/claim_status.tsv
  # attributes producers as `file.Rmd::chunk-name`, so a label makes the emission and the ledger
  # directly comparable, and unlike a line number it does not drift when the Rmd is edited above.
  src_line <- ""
  calls <- sys.calls()
  if (length(calls) > 1L) {
    src <- getSrcref(calls[[length(calls) - 1L]])
    if (!is.null(src)) src_line <- as.character(src[1L])
  }
  if (!nzchar(src_line) && requireNamespace("knitr", quietly = TRUE)) {
    lbl <- tryCatch(knitr::opts_current$get("label"), error = function(e) NULL)
    if (!is.null(lbl) && nzchar(lbl)) src_line <- paste0("chunk:", lbl)
  }

  # THE HEADER MUST MATCH tools/claim_emit.py's _HEADER EXACTLY, INCLUDING THE SIX §5 COLUMNS
  # THIS SIDE NEVER FILLS. Both emitters append to ONE file, and whichever runs first writes the
  # header. §1-4 is entirely R, so R writes it -- and with a 10-column header the Python side's
  # six appended fields land in csv.DictReader's unnamed restkey: measured, `ensemble_member`
  # then reads None, the S4 key silently degrades to (claim_id, quantity, "") and all five
  # members of a §5 quantity flag as duplicate writers. That is precisely the failure C50
  # extended the schema in place to avoid, arriving through the back door of the OTHER emitter.
  # Parity is therefore not tidiness; it is the condition under which the §5 contract holds.
  header <- c("claim_id", "quantity", "value", "published", "n", "population",
              "producer_file", "producer_line", "emitted_at", "run_id",
              "ensemble_member", "aggregation", "n_background",
              "replicate", "sd_within", "sd_between", "restriction")

  # THE CONTRACT CHECKS ITSELF, AND THE CHECK TRAVELS WITH THE FILE. This is the R half of
  # build_results_view's EMITTER_HEADER_FINGERPRINT. That guard hashes the PYTHON header and
  # never reads this file, so before this block two copies of this emitter could diverge with
  # nothing detecting it -- which was the unexamined premise in keeping a second copy at all.
  #
  # It is spelled as a literal comparison rather than a hash ON PURPOSE. A hash tells you the
  # contract moved; this tells you WHICH COLUMN, which is the difference between the failure
  # being legible and merely red. The 16-vs-17 divergence ran five days behind a test that
  # asserted `length(row) == 16` and so could not name `restriction` as the missing field.
  #
  # It lives INSIDE the function rather than beside a package import because eight Stage A steps
  # source this file inside a container that binds only <repo>:/work and <deposit>:/deposit.
  # Anything reaching outside those two mounts is unreachable there and on Explorer, where the
  # repos live under ~/cc. A check that cannot run in the environment that matters is not a check.
  expected <- c("claim_id", "quantity", "value", "published", "n", "population",
                "producer_file", "producer_line", "emitted_at", "run_id",
                "ensemble_member", "aggregation", "n_background",
                "replicate", "sd_within", "sd_between", "restriction")
  if (!identical(header, expected)) {
    stop(sprintf(paste0("claim_emit.R: this copy's header has diverged from the recorded ",
                        "contract.\n  missing: %s\n  unexpected: %s\n",
                        "The Python and R emitters write ONE store; a value emitted under a ",
                        "different header is a different artifact. Reconcile against ",
                        "claim_tools/R/claim_emit.R and bump DERIVATION_VERSION."),
                 paste(setdiff(expected, header), collapse = ", "),
                 paste(setdiff(header, expected), collapse = ", ")))
  }
  opt <- function(x) if (is.null(x)) "" else clean(x)
  row <- c(clean(claim_id), clean(quantity), clean(value),
           if (is.null(published)) "" else clean(published),
           if (is.null(n)) "" else clean(n), clean(population),
           clean(src_file), clean(src_line),
           format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
           clean(if (is.null(run_id)) Sys.getenv("NMD_RUN_ID", unset = "") else run_id),
           opt(ensemble_member), opt(aggregation), opt(n_background),
           opt(replicate), opt(sd_within), opt(sd_between),
           # NEVER BLANK, matching claim_tools/src/claim_emit.py:225. Absence is a distinct state
           # from "considered and unrestricted", and the whole point of the field is that absence
           # must not be readable as permission. An R-emitted row that left this empty reached the
           # store missing the field the fail-closed interpretability default is computed from --
           # it did not read as unrestricted, it read as absent. W1090.
           if (is.null(restriction)) "unknown" else clean(restriction))
  # A row that does not match its header corrupts every downstream read, and does it silently.
  if (length(row) != length(header))
    stop(sprintf("claim_emit: row has %d fields but the header has %d", length(row), length(header)))

  fresh <- !file.exists(target) || file.info(target)$size == 0
  con <- file(target, open = "a")
  on.exit(close(con))
  if (fresh) writeLines(paste(header, collapse = "\t"), con)
  writeLines(paste(row, collapse = "\t"), con)
  invisible(row)
}
