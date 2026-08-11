#!/usr/bin/env Rscript
# =============================================================================
# 05s_b_orfik_scan_extend.R — extend an existing ORFik scan to cover the
# isoforms it is missing, instead of rescanning everything.
#
# WHY THIS EXISTS. The model's isoform universe descends from an ORF scan run
# 2026-04-29 covering 43,039 transcripts. Measured against the CURRENT 4-CT
# labels intersected with coding status (42,063 eligible), that scan already
# holds 41,589 of them and is missing 474. A full rescan of ~42k transcripts to
# recover 474 is the single largest compute step in the §5 rerun (W32) and it
# was blocking the whole thing; extending costs minutes.
#
# WHY EXTENDING IS VALID, and this is the part that had to be checked before
# relying on it: 05s_orfik_scan.R is PER-TRANSCRIPT. Its two loops iterate one
# transcript at a time, and every aggregate inside them -- cumsum(exon_lengths),
# sum(all_atg_pos < os), sum(junctions > oe), the n_orfs / n_nmd_competent /
# n_strong_kozak counts -- is computed within a single transcript's own ORF set.
# There is no ranking, normalisation, quantile or z-score across the corpus, so
# scanning a subset and appending gives the same rows the full scan would.
# If that ever stops being true this script is WRONG, so it re-checks the scan
# parameters against the base scan's own metadata and refuses to merge on any
# mismatch.
#
# HOW IT RESTRICTS THE SCAN. It does NOT reimplement or edit 05s_orfik_scan.R --
# duplicated scan logic is how two answers to one question appear. It writes a
# FASTA containing only the missing records and re-invokes 05s with NMD_SQANTI
# pointed at it, so 05s's own `target_seqs <- all_seqs[names(all_seqs) %in%
# target_ids]` naturally reduces to that subset. Config resolution order is env
# var -> paths.yml (see REPRODUCTION.md), so no code change is needed.
#
# Usage:  Rscript analysis/isopair/isopair_wrapper/05s_b_orfik_scan_extend.R \
#             [--base <scan.rds>] [--dry-run]
#
# Writes:  .P$OUT/data_mashr/analysis_cache/orfik_scan.rds   (merged)
#          .P$OUT/data_mashr/analysis_cache/orfik_scan_extension_<n>.rds
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(Biostrings) })

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)

args     <- commandArgs(trailingOnly = TRUE)
dry_run  <- "--dry-run" %in% args
base_arg <- if ("--base" %in% args) args[which(args == "--base") + 1L] else NA_character_

CACHE     <- file.path(.P$OUT, "data_mashr/analysis_cache")
BASE_SCAN <- if (!is.na(base_arg)) base_arg else file.path(CACHE, "orfik_scan.rds")
SCAN_SRC  <- file.path(.nmd_root, "analysis/isopair/isopair_wrapper/05s_orfik_scan.R")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

stopifnot("base scan not found -- pass --base <path>" = file.exists(BASE_SCAN))
cat("=== ORFik scan extension ===\nbase scan: ", BASE_SCAN, "\n", sep = "")

# ---- 1. the eligible population, exactly as 05s defines it -------------------
# 05s:55  target_ids <- intersect(c(nmd_ids, non_nmd_ids), coding_ids)
# Recomputed here from the same inputs so the two cannot drift apart silently.
# READ THE SAME TREE 05s READS -- the requirement is unchanged, the tree is not.
#
# 05s used to setwd(.P$ISOPAIR) and read "data_mashr/..." relative to it. Computing `eligible`
# from a DIFFERENT root than 05s uses is what made the first run hand 05s 474 records and scan
# only 25: 05s's own intersection threw the rest out. Same source, or the delta is meaningless.
#
# Both now read .P$ISOPAIR_OUT (repointed 2026-07-27, Pete's decision -- see the header of
# 05s_orfik_scan.R). This must track 05s exactly; if that file's input root changes, change it
# here in the same commit. The guard below is what makes a drift loud rather than silent.
DM        <- .P$ISOPAIR_OUT
nmd_class <- readRDS(file.path(DM, "nmd_classification.rds"))
cds       <- as.data.table(readRDS(file.path(DM, "cds.rds")))
coding    <- unique(cds[coding_status == "coding", isoform_id])
eligible  <- intersect(c(nmd_class$all_samples$nmd, nmd_class$all_samples$non_nmd), coding)

base    <- readRDS(BASE_SCAN)
have    <- unique(base$tx_summary$isoform_id)
missing <- setdiff(eligible, have)

# CONVERGENCE. Without this the script could never finish: an eligible transcript with no
# ATG-initiated ORF >= 9 nt is never added to tx_summary, so it stays in `missing` forever and
# every run re-attempts the identical batch and dies the identical death. Measured 2026-07-29:
# eligible 42,063, scanned 42,043, MISSING 20 -- and all 20 yield zero ORFs, so the gap is not
# work left undone, it is a category the scan had no way to record.
#
# `no_orfs` is that category, persisted in the scan's own metadata: transcripts that WERE
# scanned and legitimately produced nothing. Subtracting them is what lets "nothing to do"
# ever fire.
known_no_orfs <- base$metadata$no_orfs
if (is.null(known_no_orfs)) known_no_orfs <- character(0)
if (length(intersect(missing, known_no_orfs))) {
  cat(sprintf("  %d of the missing were already scanned and have NO ORFs; not retrying.\n",
              length(intersect(missing, known_no_orfs))))
  missing <- setdiff(missing, known_no_orfs)
}
cat(sprintf("eligible %d | already scanned %d | recorded ORF-less %d | MISSING %d\n",
            length(eligible), length(intersect(eligible, have)),
            length(known_no_orfs), length(missing)))
if (!length(missing)) { cat("nothing to do -- base scan already covers the population\n"); quit(status = 0) }

# A transcript with no structure cannot be scanned; say so rather than let 05s drop it silently.
st <- as.data.table(readRDS(file.path(DM, "structures.rds")))
no_struct <- setdiff(missing, st$isoform_id)
if (length(no_struct)) {
  cat(sprintf("NOTE: %d missing isoform(s) have no structure and cannot be scanned; excluded.\n",
              length(no_struct)))
  missing <- setdiff(missing, no_struct)
}
if (dry_run) { cat("--dry-run: stopping before the scan\n"); quit(status = 0) }

# ---- 2. a FASTA of just those records ---------------------------------------
work <- file.path(.P$OUT, "orfik_scan_extend"); dir.create(work, recursive = TRUE, showWarnings = FALSE)
sub_dir <- file.path(work, "sqanti"); dir.create(sub_dir, showWarnings = FALSE)
sub_fa  <- file.path(sub_dir, "nmd_lungcells_corrected.fasta")
cat("extracting", length(missing), "sequences...\n")
# Streamed, not readDNAStringSet(): the full corrected FASTA is ~1.7 GB and this
# runs on an 8 GB laptop. Reading it whole to keep 474 records is how that
# machine was made unusable once already.
con_in <- file(file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta"), "r")
con_out <- file(sub_fa, "w"); want <- new.env(hash = TRUE)
for (m in missing) assign(m, TRUE, envir = want)
keep <- FALSE; n_out <- 0L
repeat {
  ln <- readLines(con_in, n = 1L, warn = FALSE); if (!length(ln)) break
  if (startsWith(ln, ">")) {
    id <- sub("^>", "", strsplit(ln, "[ \t]")[[1]][1])
    keep <- exists(id, envir = want, inherits = FALSE)
    if (keep) n_out <- n_out + 1L
  }
  if (keep) writeLines(ln, con_out)
}
close(con_in); close(con_out)
cat(sprintf("  wrote %d records -> %s\n", n_out, sub_fa))
stopifnot("not every missing isoform was found in the FASTA" = n_out == length(missing))

# ---- 3. run 05s UNMODIFIED against that subset -------------------------------
ext_out <- file.path(work, "out")
dir.create(file.path(ext_out, "data_mashr/analysis_cache"), recursive = TRUE, showWarnings = FALSE)
cat("running 05s_orfik_scan.R restricted to the subset...\n")
rc <- system2("Rscript", shQuote(SCAN_SRC),
              env = c(sprintf("NMD_SQANTI=%s", sub_dir), sprintf("NMD_OUT=%s", ext_out)),
              stdout = file.path(work, "scan.log"), stderr = file.path(work, "scan.log"))
# THE ALL-ORF-LESS BATCH IS A RESULT, NOT A FAILURE -- and telling them apart is the whole
# reason this script could not converge. 05s raises NO_ORFS_IN_BATCH when every transcript in
# the batch yields zero ORFs. That is the terminal state for these transcripts: rescanning
# cannot change it, because findORFs is deterministic on a fixed sequence and fixed parameters
# (verified independently, C45). So record them and finish, instead of failing here forever.
scan_log <- readLines(file.path(work, "scan.log"), warn = FALSE)
if (rc != 0 && any(grepl("NO_ORFS_IN_BATCH", scan_log, fixed = TRUE))) {
  cat(sprintf("\n%d transcript(s) scanned and found to have NO ORFs under start_codon=%s, ",
              length(missing), base$metadata$start_codon),
      sprintf("minimumLength=%s.\n", base$metadata$min_orf_length), sep = "")
  cat("  These are ELIGIBLE but unscannable-to-a-result. Recording them so they are not ",
      "retried, and leaving the scan's transcripts and ORFs untouched.\n", sep = "")
  merged <- base
  merged$metadata <- modifyList(base$metadata, list(
    no_orfs      = sort(unique(c(known_no_orfs, missing))),
    no_orfs_at   = Sys.time(),
    no_orfs_note = paste("scanned and returned zero ORFs; retained as a recorded category so",
                         "`missing` can reach 0 and this script can terminate")))
  if (!dry_run) saveRDS(merged, file.path(CACHE, "orfik_scan.rds"))
  cat(sprintf("  recorded %d ORF-less isoform(s); population now fully accounted for.\n",
              length(merged$metadata$no_orfs)))
  quit(status = 0)
}
if (rc != 0) { cat(tail(readLines(file.path(work, "scan.log")), 25), sep = "\n"); stop("subset scan failed") }
ext <- readRDS(file.path(ext_out, "data_mashr/analysis_cache/orfik_scan.rds"))
cat(sprintf("  scanned %d transcripts, %d ORFs\n", nrow(ext$tx_summary), nrow(ext$orf_features)))

# ---- 4. refuse to merge unless the two scans are the same experiment ---------
for (k in c("min_orf_length", "start_codon")) {
  if (!identical(base$metadata[[k]], ext$metadata[[k]]))
    stop(sprintf("parameter mismatch on %s: base=%s extension=%s -- refusing to merge",
                 k, base$metadata[[k]], ext$metadata[[k]]))
}
stopifnot("tx_summary columns differ"   = identical(names(base$tx_summary),   names(ext$tx_summary)),
          "orf_features columns differ" = identical(names(base$orf_features), names(ext$orf_features)),
          "extension overlaps the base" = !any(ext$tx_summary$isoform_id %in% base$tx_summary$isoform_id))

# ---- 4b. WHICH STRUCTURES produced each scan --------------------------------
# THE GAP THIS CLOSES (recorded 2026-07-27, unfixed until now): the checks above compare
# min_orf_length, start_codon and columns. A scan is equally a function of the exon STRUCTURES,
# and nothing recorded or compared those -- so two scans built on different structure vintages
# would have merged silently. It was benign by luck (the trees' structures proved identical
# where both have them), but the guard was weaker than the thing it guards.
#
# P5 question asked of this guard itself -- if the thing it guards against HAD happened, would
# it fire? Yes, by both routes below: a differing structures file gives a different md5, and
# differing exon structures give a different junction count per transcript.
cat("\n-- structures provenance --\n")
if (!is.null(base$metadata$structures_md5) && !is.null(ext$metadata$structures_md5)) {
  cat(sprintf("  base %s\n  ext  %s\n", base$metadata$structures_md5, ext$metadata$structures_md5))
  if (!identical(base$metadata$structures_md5, ext$metadata$structures_md5))
    stop("STRUCTURES MISMATCH: base and extension were scanned against different structures.rds",
         "\n  base: ", base$metadata$structures_path, " (", base$metadata$structures_md5, ")",
         "\n  ext : ", ext$metadata$structures_path,  " (", ext$metadata$structures_md5,  ")",
         "\n  Refusing to merge two vintages. Rescan the base against the current tree.")
  cat("  IDENTICAL -- same structures file.\n")
} else {
  # The 2026-04-29 base predates structures_md5. A hash cannot be recovered retroactively, so
  # verify by CONTENT instead, which is the thing the hash is a proxy for anyway: recompute the
  # structure-derived columns of tx_summary from the CURRENT structures and require they match
  # what the base recorded. 05s derives junctions as cumsum(exon_lengths)[-n_exons], so the count
  # is n_exons-1 (0 for single-exon); chr comes straight from structures.
  cat("  base scan predates structures_md5 (pre-2026-07-27) -- verifying BY CONTENT instead.\n")
  stopifnot("structures lacks exon_starts -- cannot verify content" = "exon_starts" %in% names(st),
            "structures lacks chr -- cannot verify content"         = "chr" %in% names(st))
  exp <- st[, .(isoform_id, chr_now = chr, njx_now = pmax(0L, lengths(exon_starts) - 1L))]
  cmp <- merge(as.data.table(base$tx_summary)[, .(isoform_id, n_junctions, chr)],
               exp, by = "isoform_id", all.x = TRUE)
  unverifiable <- cmp[is.na(njx_now)]
  checked      <- cmp[!is.na(njx_now)]
  bad_jx  <- checked[n_junctions != njx_now]
  bad_chr <- checked[as.character(chr) != as.character(chr_now)]
  cat(sprintf("  checked %s base transcripts against the current structures: %s junction mismatch, %s chr mismatch\n",
              format(nrow(checked), big.mark = ","), nrow(bad_jx), nrow(bad_chr)))
  if (nrow(bad_jx) || nrow(bad_chr)) {
    print(utils::head(rbind(bad_jx, bad_chr), 10))
    stop("STRUCTURES MISMATCH: ", nrow(bad_jx), " transcript(s) differ in junction count and ",
         nrow(bad_chr), " in chromosome between the base scan and the current structures.",
         "\n  The base was scanned against a different vintage. Refusing to merge.")
  }
  # NOT a pass for these -- they are absent from the current tree, so nothing verified them.
  cat(sprintf("  UNVERIFIABLE: %s base transcript(s) absent from the current structures.\n",
              format(nrow(unverifiable), big.mark = ",")))
}

all_tx  <- rbind(base$tx_summary,   ext$tx_summary)
all_orf <- rbind(base$orf_features, ext$orf_features)

# ---- 5. RESTRICT TO THE ELIGIBLE POPULATION ---------------------------------
# The base scan is a 6-CT-era artifact and covers transcripts the current 4-CT labels do not
# admit. Merging base+extension unfiltered carries them forward, and that is precisely how the
# mixed-provenance scan arose: 1,450 out-of-population transcripts sat in a scan in the deposit
# output tree, of which 1,365 had already reached the model's universe carrying labels the
# current classification does not assign (docs/FINDINGS.md, 2026-07-27).
#
# So the merged scan is defined as covering EXACTLY the eligible population -- the same
# intersect(labels, coding) that 05s uses. Dropping here is not a filter added to taste; it is
# what makes the scan a function of the current tree rather than of scan history.
carried <- setdiff(all_tx$isoform_id, eligible)
if (length(carried)) {
  cat(sprintf("\nDROPPING %s scanned transcript(s) not in the eligible population (6-CT-era residue).\n",
              format(length(carried), big.mark = ",")))
  cat(sprintf("  of these, %s came from the base scan and %s from the extension\n",
              format(sum(carried %in% base$tx_summary$isoform_id), big.mark = ","),
              format(sum(carried %in% ext$tx_summary$isoform_id),  big.mark = ",")))
  all_tx  <- all_tx[all_tx$isoform_id   %in% eligible, , drop = FALSE]
  all_orf <- all_orf[all_orf$isoform_id %in% eligible, , drop = FALSE]
}

merged <- list(
  orf_features = all_orf,
  tx_summary   = all_tx,
  metadata     = modifyList(base$metadata, list(
    n_transcripts   = nrow(all_tx),
    n_orfs          = nrow(all_orf),
    # provenance of the CURRENT tree, taken from the extension: the merged object is defined
    # against that tree, and a later merge must compare against it, not against the base's.
    structures_path = ext$metadata$structures_path,
    structures_md5  = ext$metadata$structures_md5,
    structures_n    = ext$metadata$structures_n,
    data_mashr_root = DM,
    restricted_to_eligible = TRUE,
    n_dropped_ineligible   = length(carried),
    extended_from   = BASE_SCAN,
    extended_n      = nrow(ext$tx_summary),
    extended_at     = Sys.time())))

saveRDS(ext,    file.path(CACHE, sprintf("orfik_scan_extension_%d.rds", nrow(ext$tx_summary))))
saveRDS(merged, file.path(CACHE, "orfik_scan.rds"))
cat(sprintf("\nMERGED: %d transcripts (%d + %d), %d ORFs\n  -> %s\n",
            nrow(merged$tx_summary), nrow(base$tx_summary), nrow(ext$tx_summary),
            nrow(merged$orf_features), file.path(CACHE, "orfik_scan.rds")))
cov <- sum(eligible %in% merged$tx_summary$isoform_id)
cat(sprintf("coverage of the eligible population: %d / %d (%.2f%%)\n",
            cov, length(eligible), 100 * cov / length(eligible)))
