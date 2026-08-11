#!/usr/bin/env Rscript
# SF26 cohort cascade — DERIVE every count the flowchart prints.
#
# Why this file exists. Until 2026-08-09 the flowchart's upstream nodes read
# nrow(profiles_c2) live while every subset box below carried a hardcoded string. The
# 2026-08-08 rebuild therefore moved half the figure onto the deposit and left the rest
# behind, and the rendered .dot ended up printing 1,548 TWICE and 1,578 ONCE -- with the
# 1,578 arrow pointing straight into the box labelled 1,548. The same box stated the
# reference-share floor yields 5,001 three lines below reporting 5,049 candidates.
#
# The literals were literals because nobody expected the substrate to move. Retyping them
# would rebuild exactly the same trap, so every number is derived here instead and
# build_flowchart.R consumes this table. A number that is not in the table is not in the
# figure.
#
# The filter chain replicates analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd
# chunks sec1-pop-bc, sec2-scope, sec2a-ptc and sec3-scope. mechanism_class() comes from
# figures/lib/mechanism_class.R rather than being inlined, so the two cannot drift.
#
# Output: data/cascade_counts.tsv  (key, value, note)

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({ library(data.table); library(Isopair) })

DM   <- nmd_data_dir(.P)
HERE <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- a[grep("--file=", a)]
  if (length(m) > 0) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  normalizePath(getwd())
})()

source(file.path(.nmd_root, "figures", "lib", "mechanism_class.R"))

expr_mat    <- readRDS(file.path(DM, "expression_data.rds"))
nmd_class   <- readRDS(file.path(DM, "nmd_classification.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
cds         <- as.data.table(readRDS(file.path(DM, "cds.rds")))
structures  <- as.data.table(readRDS(file.path(DM, "structures.rds")))
ref_atg     <- readRDS(file.path(DM, "analysis_cache/ref_atg_analysis.rds"))

out <- list()
put <- function(k, v, note = "") out[[length(out) + 1L]] <<-
  data.frame(key = k, value = v, note = note, stringsAsFactors = FALSE)

is_enst  <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

# ── Upstream ────────────────────────────────────────────────────────────────────────
put("n_raw_iso",     645273L, "raw isocall matrix rows -- a constant of the input data")
put("n_filter_iso",  nrow(expr_mat))
put("n_4ct_samples", 26L, "AT 6 + DD 8 + FB 6 + MV 6")
put("n_nmd_as",      length(nmd_class[["all_samples"]]$nmd))
put("n_nonnmd_as",   length(nmd_class[["all_samples"]]$non_nmd))
put("n_c2",          nrow(profiles_c2))
put("n_c4",          nrow(profiles_c4))

# The reference-share floor's PRE-floor count. 02_build_profiles_mashr.R now EMITS it as data --
# reference_floor_counts_allsamples.tsv beside the caches -- because it exists nowhere else: both
# pairs_c4 and profiles_c4 are written post-floor, so the number the flowchart prints was only ever
# in stdout.
#
# THIS READ A LOG UNTIL 2026-08-09 AND THE CLEAN ROOM KILLED IT. Job 9045389 on Explorer stopped
# here with "missing /work/tmp/out/logs/02_floored.log". A log is not an artifact: it is not carried
# between jobs, not deposited, and not reproducible by anyone who runs step 02 and this figure in
# separate sittings -- which is exactly what a reader does. The failure was invisible to every
# static check because the file existed on the authoring laptop.
#
# The log path is kept as a FALLBACK for trees built before the emitter existed, and only as that.
# Either way the post-floor count is cross-checked against nrow(profiles_c4), so a stale source
# stops the build instead of printing a mismatched pair.
floor_tsv <- file.path(nmd_data_dir(.P), "reference_floor_counts_allsamples.tsv")
floor_log <- file.path(.P$OUT, "logs", "02_floored.log")
if (file.exists(floor_tsv)) {
  ft <- utils::read.delim(floor_tsv, stringsAsFactors = FALSE)
  gv <- function(k) as.integer(ft$value[match(k, ft$metric)])
  floor_before <- gv("c4_pairs_before_floor"); floor_after <- gv("c4_pairs_after_floor")
  floor_src <- "reference_floor_counts_allsamples.tsv"
} else if (file.exists(floor_log)) {
  fl <- grep("Reference-share floor", readLines(floor_log, warn = FALSE), value = TRUE)[1]
  if (is.na(fl)) stop("no reference-share floor line in ", floor_log)
  m <- regmatches(fl, regexec("([0-9]+) -> ([0-9]+) C4", fl))[[1]]
  if (length(m) != 3L) stop("cannot parse the floor line: ", fl)
  floor_before <- as.integer(m[2]); floor_after <- as.integer(m[3])
  floor_src <- "02_floored.log (LEGACY fallback -- rerun 02 to emit the tsv)"
} else {
  stop("no reference-floor counts: expected ", floor_tsv,
       "\n  Rerun analysis/isopair/isopair_wrapper/02_build_profiles_mashr.R, which emits it.")
}
floor_dropped <- floor_before - floor_after
if (floor_after != nrow(profiles_c4))
  stop(sprintf(paste("%s reports %d C4 pairs after the floor but profiles_c4_allsamples.rds holds",
                     "%d -- different vintages. Rerun 02 before rebuilding this figure."),
               floor_src, floor_after, nrow(profiles_c4)))
put("floor_before",  floor_before, "genes before the 25% reference-share floor")
put("floor_after",   floor_after,  paste0("cross-checked against nrow(profiles_c4); source: ", floor_src))
put("floor_dropped", floor_dropped)

# ── pop_BC: gene-matched on shared (gene, reference) ────────────────────────────────
pop_keys  <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_bc_c2 <- profiles_c2[make_key(profiles_c2) %in% pop_keys]
pop_bc_c4 <- profiles_c4[make_key(profiles_c4) %in% pop_keys]
put("pop_bc_c2", nrow(pop_bc_c2), "the matched-pair hub")
put("pop_bc_c4", nrow(pop_bc_c4))

# ── Section A: all-ENST -> curated CDS -> re-intersect ──────────────────────────────
c2_E <- pop_bc_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
c4_E <- pop_bc_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
j1   <- intersect(make_key(c2_E), make_key(c4_E))
c2_E <- c2_E[make_key(c2_E) %in% j1]; c4_E <- c4_E[make_key(c4_E) %in% j1]
put("secA_enst_kept",    nrow(c2_E))
put("secA_enst_dropped", nrow(pop_bc_c2) - nrow(c2_E))

coding_ids <- cds[coding_status == "coding", isoform_id]
have_cds <- function(d) d$reference_isoform_id %in% coding_ids &
                        d$comparator_isoform_id %in% coding_ids
c2_E2 <- c2_E[have_cds(c2_E)]; c4_E2 <- c4_E[have_cds(c4_E)]
j2 <- intersect(make_key(c2_E2), make_key(c4_E2))
secA_c2 <- c2_E2[make_key(c2_E2) %in% j2]; secA_c4 <- c4_E2[make_key(c4_E2) %in% j2]
put("secA_coding_kept",    nrow(secA_c2))
put("secA_coding_dropped", nrow(c2_E) - nrow(secA_c2))
put("secA_n",              nrow(secA_c2), "GENCODE-restricted subset")
stopifnot(nrow(secA_c2) == nrow(secA_c4))

# Section A PTC call: the comparator's OWN GENCODE stop, 50-nt rule.
coding_cds <- cds[coding_status == "coding"]
si <- match(coding_cds$isoform_id, structures$isoform_id)
stop_tx <- rep(NA_integer_, nrow(coding_cds))
for (i in seq_len(nrow(coding_cds))) {
  k <- si[i]; if (is.na(k)) next
  st <- coding_cds$strand[i]; if (is.na(st) || st == "") next
  g <- if (st == "+") coding_cds$cds_stop[i] else coding_cds$cds_start[i]
  stop_tx[i] <- Isopair::genomicToTranscript(
    g, structures$exon_starts[[k]], structures$exon_ends[[k]], st)
}
own_stop <- data.table(isoform_id = coding_cds$isoform_id, own_stop_tx_pos = stop_tx)

last_ejc <- function(iso) {
  k <- match(iso, structures$isoform_id); if (is.na(k)) return(NA_integer_)
  es <- structures$exon_starts[[k]]; ee <- structures$exon_ends[[k]]
  st <- structures$strand[k]; if (is.na(st) || st == "") return(NA_integer_)
  L <- ee - es + 1L; if (st == "-") L <- rev(L)
  if (length(L) <= 1L) return(NA_integer_)
  max(cumsum(L[-length(L)]))
}
ptc_call <- function(pop) {
  d <- merge(pop[, .(comparator_isoform_id)], own_stop,
             by.x = "comparator_isoform_id", by.y = "isoform_id", all.x = TRUE)
  d[, last_ejc_tx_pos := vapply(comparator_isoform_id, last_ejc, integer(1))]
  d[, has_ptc := !is.na(last_ejc_tx_pos - own_stop_tx_pos) &
                 (last_ejc_tx_pos - own_stop_tx_pos) > 50L]
  d
}
tA2 <- ptc_call(secA_c2)
put("secA_ptc_pos", sum(tA2$has_ptc, na.rm = TRUE))
put("secA_ptc_neg", nrow(secA_c2) - sum(tA2$has_ptc, na.rm = TRUE))
put("secA_ptc_pct", round(100 * sum(tA2$has_ptc, na.rm = TRUE) / nrow(secA_c2)))

# ── Section C: ENST reference -> ref-AUG-traceable -> re-intersect ──────────────────
c2_base <- pop_bc_c2[is_enst(reference_isoform_id)]
c4_base <- pop_bc_c4[is_enst(reference_isoform_id)]
put("secC_enstref_kept",    nrow(c2_base))
put("secC_enstref_dropped", nrow(pop_bc_c2) - nrow(c2_base))

ra2 <- as.data.table(ref_atg$c2); ra4 <- as.data.table(ref_atg$c4)
c2_full <- ra2[comparator_isoform_id %in% c2_base$comparator_isoform_id]
c4_full <- ra4[comparator_isoform_id %in% c4_base$comparator_isoform_id]
c2_tr <- c2_full[category %in% REF_AUG_PRESENT_CATEGORIES]
c4_tr <- c4_full[category %in% REF_AUG_PRESENT_CATEGORIES]
put("secC_traceable_kept_nmd",     nrow(c2_tr))
put("secC_traceable_kept_ctrl",    nrow(c4_tr))
# Dropped is measured against the PREVIOUS BOX (the ENST-reference set), not against
# c2_full -- ref_atg does not cover every ENST-reference comparator, and counting only the
# category filter would leave kept + dropped short of the box above it. A flowchart whose
# rows do not sum to the row above is unreadable, so the arithmetic is asserted below.
put("secC_traceable_dropped_nmd",  nrow(c2_base) - nrow(c2_tr))
put("secC_traceable_dropped_ctrl", nrow(c4_base) - nrow(c4_tr))

sk <- merge(unique(c2_tr[, .(gene_id, reference_isoform_id)]),
            unique(c4_tr[, .(gene_id, reference_isoform_id)]),
            by = c("gene_id", "reference_isoform_id"))
secC_c2 <- merge(c2_tr, sk, by = c("gene_id", "reference_isoform_id"))
secC_c4 <- merge(c4_tr, sk, by = c("gene_id", "reference_isoform_id"))
put("secC_n",                nrow(secC_c2), "reference-AUG-traceable subset")
put("secC_match_dropped_nmd",  nrow(c2_tr) - nrow(secC_c2))
put("secC_match_dropped_ctrl", nrow(c4_tr) - nrow(secC_c4))
stopifnot(nrow(secC_c2) == nrow(secC_c4))

secC_c2[, mech4 := mechanism_class_4(
  mechanism_class(category, comp_orf_length, ref_orf_length))]
n_pp <- sum(secC_c2$mech4 == "NMD+/PTC+")
put("secC_ptc_pos", n_pp)
put("secC_ptc_neg", nrow(secC_c2) - n_pp)
put("secC_ptc_pct", round(100 * n_pp / nrow(secC_c2)))

# ── Occult-PTC: ref-AUG says premature stop, TD2's own CDS call does not ────────────
occ <- secC_c2[category == "effectively_ptc" & original_ptc == FALSE]
put("occult_n",       nrow(occ), "occult-PTC subset")
put("occult_dropped", nrow(secC_c2) - nrow(occ), "no disagreement between the two CDS calls")

res <- do.call(rbind, out)

# ── Every cascade box must sum to the box above it ──────────────────────────────────
# The defect this figure carried was two boxes disagreeing about one population, so the
# rebuild asserts the arithmetic rather than trusting it. Each check names the two boxes
# a reader would compare by eye.
g <- function(k) res$value[match(k, res$key)]
chk <- function(label, lhs, rhs) if (lhs != rhs)
  stop(sprintf("%s: %d != %d -- the flowchart would print boxes that do not add up",
               label, lhs, rhs))
chk("pop_BC arms",        g("pop_bc_c2"),       g("pop_bc_c4"))
chk("hub == C2",          g("pop_bc_c2"),       g("n_c2"))
chk("floor arithmetic",   g("floor_before") - g("floor_dropped"), g("floor_after"))
chk("secA ENST step",     g("secA_enst_kept") + g("secA_enst_dropped"),       g("pop_bc_c2"))
chk("secA coding step",   g("secA_coding_kept") + g("secA_coding_dropped"),   g("secA_enst_kept"))
chk("secA groups",        g("secA_ptc_pos") + g("secA_ptc_neg"),              g("secA_n"))
chk("secC ENST-ref step", g("secC_enstref_kept") + g("secC_enstref_dropped"), g("pop_bc_c2"))
chk("secC traceable NMD", g("secC_traceable_kept_nmd") + g("secC_traceable_dropped_nmd"),
                          g("secC_enstref_kept"))
chk("secC traceable ctrl", g("secC_traceable_kept_ctrl") + g("secC_traceable_dropped_ctrl"),
                           g("secC_enstref_kept"))
chk("secC match NMD",     g("secC_n") + g("secC_match_dropped_nmd"),  g("secC_traceable_kept_nmd"))
chk("secC match ctrl",    g("secC_n") + g("secC_match_dropped_ctrl"), g("secC_traceable_kept_ctrl"))
chk("secC groups",        g("secC_ptc_pos") + g("secC_ptc_neg"),      g("secC_n"))
chk("occult step",        g("occult_n") + g("occult_dropped"),        g("secC_n"))
cat("[check] all 13 cascade sums close\n")

dir.create(file.path(HERE, "data"), showWarnings = FALSE)
tsv <- file.path(HERE, "data", "cascade_counts.tsv")
write.table(res, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[tsv]  -> %s  (%d counts)\n", tsv, nrow(res)))
print(res, row.names = FALSE)
