# =============================================================================
# mashr cell-type column resolution — SSOT
#
# The shared mashr objects exist under two cell-type conventions:
#   internal codes  Smg1i_in_AT   Smg1i_in_DD    (pre April-2026 rename)
#   published names Smg1i_in_AT2  Smg1i_in_LAE   (current; what Yul's code expects)
#
# FB and MV are spelled identically in both, so a mismatch resolves only two of four
# cell types. The failure mode is therefore a SILENT ZERO, not an error -- on 2026-07-25
# it dropped a reported range from 20.7-54.5% to 0-38.9% without any warning.
#
# Resolve through here rather than hard-coding either spelling.
# =============================================================================

# Candidate column names for one cell type, preferred spelling first.
nmd_mash_candidates <- function(ct) {
  switch(toupper(ct),
    AT2 = c("Smg1i_in_AT2", "Smg1i_in_AT"),
    AT  = c("Smg1i_in_AT2", "Smg1i_in_AT"),
    LAE = c("Smg1i_in_LAE", "Smg1i_in_DD"),
    DD  = c("Smg1i_in_LAE", "Smg1i_in_DD"),
    FB  = "Smg1i_in_FB",
    MV  = "Smg1i_in_MV",
    paste0("Smg1i_in_", ct))
}

# Resolve one cell type against one or more data frames. Stops loudly if absent.
nmd_mash_col <- function(ct, ...) {
  frames <- list(...)
  stopifnot("supply at least one data frame" = length(frames) > 0)
  cand <- nmd_mash_candidates(ct)
  # colnames(), not names(): ashr::get_pm() and get_lfsr() return MATRICES, and
  # names(matrix) is NULL -- so every candidate missed and the resolver stopped with
  # "no mashr column". colnames() is correct for both matrices and data frames.
  ok <- Reduce(intersect, lapply(frames, colnames))
  hit <- cand[cand %in% ok]
  if (!length(hit)) {
    stop("no mashr column for cell type '", ct, "' -- tried: ", paste(cand, collapse = ", "),
         "\n  available: ", paste(grep("^Smg1i_in_", ok, value = TRUE), collapse = ", "))
  }
  hit[1]
}

# Resolve the four manuscript cell types at once, in the given order.
nmd_mash_cols <- function(cts = c("AT2", "LAE", "FB", "MV"), ...) {
  vapply(cts, nmd_mash_col, character(1), ...)
}

# Rename whatever spelling a frame carries to the canonical PUBLISHED names, so downstream
# code can use Smg1i_in_AT2 / Smg1i_in_LAE unconditionally. Normalising once at load is
# preferable to resolving at every use site.
nmd_normalize_mash_cols <- function(df) {
  nm <- names(df)
  nm[nm == "Smg1i_in_AT"] <- "Smg1i_in_AT2"
  nm[nm == "Smg1i_in_DD"] <- "Smg1i_in_LAE"
  names(df) <- nm
  df
}

# =============================================================================
# FILENAME resolution — the same two-convention problem, one layer down.
#
# The per-cell-type mashr exports carry the SAME two spellings as the columns:
#   post-rename   nmd_mashr_die_at2_*.csv   nmd_mashr_die_lae_*.csv
#   pre-rename    nmd_mashr_die_at_*.csv    nmd_mashr_die_dd_*.csv
# FB and MV are identical in both, so a mismatch resolves EXACTLY TWO OF FOUR -- and the
# call sites did `if (!file.exists(f)) return(NULL)`, so the other two vanished in silence.
# On 2026-07-26 build_paper_tables.R reported "2 cell type rows" for a four-cell-type table
# and nothing anywhere said a word.
#
# That is the identical failure this file was written to stop, one layer down from the
# columns it already guarded. Resolve through here and it STOPS instead of half-answering.
# =============================================================================
nmd_mash_file_candidates <- function(ct) {
  switch(toupper(ct),
    AT2 = c("at2", "at"), AT = c("at2", "at"),
    LAE = c("lae", "dd"),  DD = c("lae", "dd"),
    FB  = "fb", MV = "mv",
    tolower(ct))
}

# Resolve one per-cell-type export. `pattern` takes the suffix, e.g.
#   nmd_mash_file(lr_dir, function(s) sprintf("nmd_mashr_die_%s_2026.3.10.csv", s), "DD")
nmd_mash_file <- function(dir, pattern, ct, required = TRUE) {
  cand <- nmd_mash_file_candidates(ct)
  for (s in cand) {
    f <- file.path(dir, pattern(s))
    if (file.exists(f)) return(f)
  }
  if (!required) return(NA_character_)
  stop("no mashr export for cell type '", ct, "' in\n    ", dir,
       "\n  tried: ", paste(vapply(cand, pattern, character(1)), collapse = ", "),
       "\n  A missing file here silently drops one cell type from a four-cell-type table.",
       call. = FALSE)
}
