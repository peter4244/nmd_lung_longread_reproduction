#!/usr/bin/env Rscript
# Consolidated TD2-bias evidence supplement.
#
# Scope (matches Figure 4 Panels C/D's universe — 1:1 gene-matched, re-intersected):
#   pop_BC NMD c2 + Control c4 sharing (gene_id, reference_isoform_id)
#   ∩ ENST-reference (reference_isoform_id starts with "ENST")
#   ∩ ref-AUG-traceable category (effectively_ptc / no_downstream_ejc / truncated_no_ejc) on EACH side
#   ∩ re-intersect (gene_id, reference_isoform_id) AFTER the category filter
#   ∩ TD2 ORF length computable (cds.rds orf_length not NA)
#   ∩ ref-AUG-anchored ORF length computable
#   ∩ ref AUG exonic in comparator
#   → n = 1,166 NMD comparator pairs (Row 1)
#
# Row 2 occult-PTC subset (effectively_ptc ∩ original_ptc == FALSE within the
# re-intersected NMD universe):
#   → n = 492 NMD comparator pairs
#
# Within this universe, ~53% of pairs have TD2 selecting a different ATG than
# the reference AUG (the population in which TD2's start-codon choice differs
# from the reference and can therefore be compared). The remaining ~47% have
# TD2 choosing the same start codon as the reference AUG, in which case all
# three quantities (TD2 ORF length, Kozak score, position) trivially match
# the reference side — these pairs add no signal but do anchor the denominator.
#
# Three observations on this set:
#   Panel A — TD2-selected CDS ORF length vs reference-AUG-anchored ORF length
#   Panel B — Kozak PWM at the reference AUG vs at the TD2-predicted ATG (paired)
#   Panel C — TD2 ATG position relative to the reference AUG (transcript coords)

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(data.table)
})

# Repointed 2026-07-26 (W31 step 2). Was .P$ISOPAIR -- the LEGACY tree -- so this

# exported panels from the 2026-03-11 structure universe while the reports use the

# rebuilt one. Resolver lives in R/load_config.R; --data-dir still overrides.

ISOPAIR_DIR <- dirname(nmd_data_dir(.P))
DM <- file.path(ISOPAIR_DIR, "data_mashr")
FASTA_PATH <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")

ref_atg <- readRDS(file.path(DM, "analysis_cache/ref_atg_analysis.rds"))
cds     <- as.data.table(readRDS(file.path(DM, "cds.rds")))
struct  <- as.data.table(readRDS(file.path(DM, "structures.rds")))
prof_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
prof_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

c2 <- as.data.table(ref_atg$c2)
c4 <- as.data.table(ref_atg$c4)
cat(sprintf("ref_atg$c2 universe: n = %d\n", nrow(c2)))

# Build the re-intersected NMD universe (matches Fig 4 Section C scope).
# Step 1: pop_BC (1:1 gene-matched at the c2 + c4 level).
keys_BC <- merge(
  unique(prof_c2[, .(gene_id, reference_isoform_id)]),
  unique(prof_c4[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
# Step 2: restrict ref_atg c2 + c4 to ENST-reference + pop_BC + ref-AUG-traceable.
is_enst_fn <- function(x) startsWith(x, "ENST")
TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
c2_tr <- merge(c2[is_enst_fn(reference_isoform_id) & category %in% TR],
               keys_BC, by = c("gene_id", "reference_isoform_id"))
c4_tr <- merge(c4[is_enst_fn(reference_isoform_id) & category %in% TR],
               keys_BC, by = c("gene_id", "reference_isoform_id"))
# Step 3: re-intersect on (gene_id, reference_isoform_id) AFTER category filter.
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
c2 <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
cat(sprintf("After ENST-ref + TR + re-intersect with Control side: n = %d\n",
            nrow(c2)))

# Add comparator CDS metadata
c2_cds <- cds[match(c2$comparator_isoform_id, isoform_id)]
c2[, strand_use := c2_cds$strand]
c2[, td2_atg_genomic := ifelse(c2_cds$strand == "+",
                                c2_cds$cds_start, c2_cds$cds_stop)]
c2[, td2_orf_nt := c2_cds$orf_length * 3L]
c2[, ref_atg_orf_nt := comp_orf_length_from_ref_atg]

# Additional filters: TD2 CDS + ref-AUG ORF + ref AUG exonic (so all 3 panels
# can be measured on the same isoforms).
broad <- c2[!is.na(td2_orf_nt) & !is.na(ref_atg_orf_nt) &
            ref_atg_exonic_in_comp == TRUE]
cat(sprintf("After ENST-ref + ref-AUG-traceable + TD2 CDS + ref-AUG ORF + ref AUG exonic: n = %d\n",
            nrow(broad)))
cat(sprintf("  Of which TD2 ATG != ref AUG: %d (%.1f%%)\n",
            sum(broad$td2_atg_genomic != broad$ref_atg_genomic, na.rm = TRUE),
            100 * mean(broad$td2_atg_genomic != broad$ref_atg_genomic, na.rm = TRUE)))

# Genomic -> comparator transcript coord converter
genomic_to_tx <- function(gp, es, ee, s) {
  if (s == "-") { es <- rev(es); ee <- rev(ee) }
  cum <- 0L
  for (j in seq_along(es)) {
    if (s == "+" && gp >= es[j] && gp <= ee[j]) return(cum + (gp - es[j]) + 1L)
    if (s == "-" && gp >= es[j] && gp <= ee[j]) return(cum + (ee[j] - gp) + 1L)
    cum <- cum + (ee[j] - es[j] + 1L)
  }
  NA_integer_
}
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

broad_struct <- struct[match(broad$comparator_isoform_id, isoform_id)]
ref_tx <- integer(nrow(broad))
td2_tx <- integer(nrow(broad))
for (i in seq_len(nrow(broad))) {
  es <- broad_struct$exon_starts[[i]]; ee <- broad_struct$exon_ends[[i]]
  s  <- broad$strand_use[i]
  if (length(es) == 0L || is.na(broad$ref_atg_genomic[i]) ||
      is.na(broad$td2_atg_genomic[i])) {
    ref_tx[i] <- NA_integer_; td2_tx[i] <- NA_integer_; next
  }
  ref_tx[i] <- genomic_to_tx(broad$ref_atg_genomic[i], es, ee, s) %||% NA_integer_
  td2_tx[i] <- genomic_to_tx(broad$td2_atg_genomic[i], es, ee, s) %||% NA_integer_
}
broad[, ref_atg_tx_pos := ref_tx]
broad[, td2_atg_tx_pos := td2_tx]
broad[, distance_nt   := td2_tx - ref_tx]

# Keep only pairs with both ATGs mappable
keep <- !is.na(broad$ref_atg_tx_pos) & !is.na(broad$td2_atg_tx_pos)
broad_clean <- broad[keep]
cat(sprintf("After both ATGs mappable to comp tx coords: n = %d\n", nrow(broad_clean)))

# Repointed 2026-07-26 (W31 step 2): was the LEGACY figures tree (figures/SupplementalFigures/, renamed `supplemental/` here).
SUPP_DIR <- file.path(.P$ROOT, "figures", "supplemental", "SF34_TD2Bias_broad")
dir.create(file.path(SUPP_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Panel A: TD2 vs reference-AUG ORF length
# ─────────────────────────────────────────────────────────────────────────────
panelA <- data.table(
  comparator_isoform_id = rep(broad_clean$comparator_isoform_id, 2),
  ORF_source = rep(c("TD2-selected CDS", "Reference ATG ORF"),
                   each = nrow(broad_clean)),
  length_nt  = c(broad_clean$td2_orf_nt, broad_clean$ref_atg_orf_nt))
fwrite(panelA, file.path(SUPP_DIR, "data", "panelA_td2_vs_refaug_length.tsv"),
       sep = "\t")
cat(sprintf("\n[Panel A] TD2 vs ref-AUG ORF length: n=%d each\n", nrow(broad_clean)))
cat(sprintf("  TD2 median = %.0f, ref-AUG median = %.0f, ratio = %.2fx\n",
            median(broad_clean$td2_orf_nt), median(broad_clean$ref_atg_orf_nt),
            median(broad_clean$td2_orf_nt) / median(broad_clean$ref_atg_orf_nt)))

# ─────────────────────────────────────────────────────────────────────────────
# Panel B: Paired Kozak PWM scores at reference AUG vs TD2-predicted ATG
# ─────────────────────────────────────────────────────────────────────────────
load_sequences <- function(isoform_ids) {
  ids <- unique(isoform_ids)
  id_file <- tempfile(fileext = ".txt")
  writeLines(ids, id_file)
  on.exit(unlink(id_file), add = TRUE)
  awk_cmd <- sprintf(
    "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
    id_file, FASTA_PATH)
  raw <- system(awk_cmd, intern = TRUE)
  split_raw <- strsplit(raw, "\t", fixed = TRUE)
  setNames(toupper(vapply(split_raw, `[`, character(1), 2)),
                   vapply(split_raw, `[`, character(1), 1))
}
seqs <- load_sequences(broad_clean$comparator_isoform_id)
pwm_seqs <- character(nrow(broad_clean))
for (k in seq_len(nrow(broad_clean))) {
  cid <- broad_clean$comparator_isoform_id[k]
  if (cid %in% names(seqs)) pwm_seqs[k] <- seqs[cid]
}
v_ref <- !is.na(broad_clean$ref_atg_tx_pos) & pwm_seqs != ""
v_td2 <- !is.na(broad_clean$td2_atg_tx_pos) & pwm_seqs != ""
ref_pwm <- rep(NA_real_, nrow(broad_clean))
td2_pwm <- rep(NA_real_, nrow(broad_clean))
if (sum(v_ref) > 0) {
  res <- Isopair::scoreKozakPWM(pwm_seqs[v_ref], broad_clean$ref_atg_tx_pos[v_ref])
  ref_pwm[v_ref] <- res$score
}
if (sum(v_td2) > 0) {
  res <- Isopair::scoreKozakPWM(pwm_seqs[v_td2], broad_clean$td2_atg_tx_pos[v_td2])
  td2_pwm[v_td2] <- res$score
}
koz_valid <- !is.na(ref_pwm) & !is.na(td2_pwm)
n_kozak <- sum(koz_valid)
wt <- wilcox.test(ref_pwm[koz_valid], td2_pwm[koz_valid], paired = TRUE)
panelB <- data.table(
  comparator_isoform_id = rep(broad_clean$comparator_isoform_id[koz_valid], 2),
  ATG_source = rep(c("Reference ATG", "TD2-predicted ATG"),
                   each = n_kozak),
  score = c(ref_pwm[koz_valid], td2_pwm[koz_valid]))
fwrite(panelB, file.path(SUPP_DIR, "data", "panelB_kozak.tsv"), sep = "\t")
cat(sprintf("\n[Panel B] Kozak PWM: n=%d  ref median=%.2f  TD2 median=%.2f  paired-Wilcoxon p=%.3g\n",
            n_kozak, median(ref_pwm[koz_valid]),
            median(td2_pwm[koz_valid]), wt$p.value))

# ─────────────────────────────────────────────────────────────────────────────
# Panel C: TD2 ATG position relative to reference AUG
# ─────────────────────────────────────────────────────────────────────────────
panelC <- broad_clean[, .(comparator_isoform_id, gene_id, strand = strand_use,
                          ref_atg_tx_pos, td2_atg_tx_pos, distance_nt,
                          n_downstream_ejc)]
fwrite(panelC, file.path(SUPP_DIR, "data", "panelC_td2_position.tsv"), sep = "\t")
n_eq   <- sum(panelC$distance_nt == 0L)
n_down <- sum(panelC$distance_nt > 0L)
n_up   <- sum(panelC$distance_nt < 0L)
cat(sprintf("\n[Panel C] TD2 ATG position relative to ref AUG: n=%d\n", nrow(panelC)))
cat(sprintf("  TD2 == ref AUG (distance 0):   %d (%.1f%%)\n",
            n_eq, 100 * n_eq / nrow(panelC)))
cat(sprintf("  TD2 downstream of ref AUG:     %d (%.1f%%)  median offset = +%d nt\n",
            n_down, 100 * n_down / nrow(panelC),
            if (n_down > 0) as.integer(round(median(panelC$distance_nt[panelC$distance_nt > 0]))) else 0L))
cat(sprintf("  TD2 upstream of ref AUG:       %d (%.1f%%)\n",
            n_up, 100 * n_up / nrow(panelC)))

# ─────────────────────────────────────────────────────────────────────────────
# Row 2 of the supplement: same three quantities restricted to the n=900
# "occult-PTC" subset (effectively_ptc ∩ original_ptc == FALSE — pairs in which
# reference-AUG tracing revealed a PTC that TD2 missed). This is the population
# in which TD2's PTC-avoidance bias is most pronounced.
# ─────────────────────────────────────────────────────────────────────────────
occ_clean <- broad_clean[category == "effectively_ptc" & original_ptc == FALSE]
cat(sprintf("\nOccult-PTC subset (effectively_ptc ∩ original_ptc=FALSE): n = %d\n",
            nrow(occ_clean)))

# Panel D (= Panel A on the occult-PTC subset)
panelD <- data.table(
  comparator_isoform_id = rep(occ_clean$comparator_isoform_id, 2),
  ORF_source = rep(c("TD2-selected CDS", "Reference ATG ORF"),
                   each = nrow(occ_clean)),
  length_nt  = c(occ_clean$td2_orf_nt, occ_clean$ref_atg_orf_nt))
fwrite(panelD, file.path(SUPP_DIR, "data", "panelD_td2_vs_refaug_length_occult.tsv"),
       sep = "\t")
cat(sprintf("[Panel D] (occult-PTC) TD2 vs ref-AUG ORF length: n=%d each\n",
            nrow(occ_clean)))
cat(sprintf("  TD2 median = %.0f, ref-AUG median = %.0f, ratio = %.2fx\n",
            median(occ_clean$td2_orf_nt), median(occ_clean$ref_atg_orf_nt),
            median(occ_clean$td2_orf_nt) / median(occ_clean$ref_atg_orf_nt)))

# Panel E (= Panel B on the occult-PTC subset)
occ_ids <- occ_clean$comparator_isoform_id
panelE_idx <- match(occ_ids, broad_clean$comparator_isoform_id)
ref_pwm_occ <- ref_pwm[panelE_idx]
td2_pwm_occ <- td2_pwm[panelE_idx]
koz_valid_occ <- !is.na(ref_pwm_occ) & !is.na(td2_pwm_occ)
n_kozak_occ <- sum(koz_valid_occ)
wt_occ <- wilcox.test(ref_pwm_occ[koz_valid_occ], td2_pwm_occ[koz_valid_occ],
                     paired = TRUE)
panelE <- data.table(
  comparator_isoform_id = rep(occ_ids[koz_valid_occ], 2),
  ATG_source = rep(c("Reference ATG", "TD2-predicted ATG"), each = n_kozak_occ),
  score = c(ref_pwm_occ[koz_valid_occ], td2_pwm_occ[koz_valid_occ]))
fwrite(panelE, file.path(SUPP_DIR, "data", "panelE_kozak_occult.tsv"), sep = "\t")
cat(sprintf("[Panel E] (occult-PTC) Kozak PWM: n=%d  ref median=%.2f  TD2 median=%.2f  paired-Wilcoxon p=%.3g\n",
            n_kozak_occ, median(ref_pwm_occ[koz_valid_occ]),
            median(td2_pwm_occ[koz_valid_occ]), wt_occ$p.value))

# Panel F (= Panel C on the occult-PTC subset)
panelF <- occ_clean[, .(comparator_isoform_id, gene_id, strand = strand_use,
                        ref_atg_tx_pos, td2_atg_tx_pos, distance_nt,
                        n_downstream_ejc)]
fwrite(panelF, file.path(SUPP_DIR, "data", "panelF_td2_position_occult.tsv"),
       sep = "\t")
n_eq_o   <- sum(panelF$distance_nt == 0L)
n_down_o <- sum(panelF$distance_nt > 0L)
n_up_o   <- sum(panelF$distance_nt < 0L)
cat(sprintf("[Panel F] (occult-PTC) TD2 position relative to ref AUG: n=%d\n", nrow(panelF)))
cat(sprintf("  TD2 == ref AUG: %d (%.1f%%)  TD2 downstream: %d (%.1f%%)  TD2 upstream: %d (%.1f%%)\n",
            n_eq_o, 100 * n_eq_o / nrow(panelF),
            n_down_o, 100 * n_down_o / nrow(panelF),
            n_up_o, 100 * n_up_o / nrow(panelF)))

# Summary TSV with key numbers used by the figure
summary_dt <- data.table(
  panel = c("A", "B", "C", "D", "E", "F"),
  scope = c("broad (n=1,166)", "broad (n=1,166)", "broad (n=1,166)",
            "occult-PTC (n=492)", "occult-PTC (n=492)", "occult-PTC (n=492)"),
  n_pairs = c(nrow(broad_clean), n_kozak, nrow(panelC),
              nrow(occ_clean), n_kozak_occ, nrow(panelF)),
  metric = c("TD2 vs ref-AUG ORF length",
             "Paired Kozak PWM",
             "TD2 ATG position relative to ref AUG",
             "TD2 vs ref-AUG ORF length (occult-PTC)",
             "Paired Kozak PWM (occult-PTC)",
             "TD2 ATG position relative to ref AUG (occult-PTC)"),
  headline = c(sprintf("TD2 median %.0f nt vs ref-AUG median %.0f nt (%.1fx)",
                       median(broad_clean$td2_orf_nt),
                       median(broad_clean$ref_atg_orf_nt),
                       median(broad_clean$td2_orf_nt) / median(broad_clean$ref_atg_orf_nt)),
               sprintf("ref median %.2f vs TD2 median %.2f, paired Wilcoxon p=%.3g",
                       median(ref_pwm[koz_valid]),
                       median(td2_pwm[koz_valid]), wt$p.value),
               sprintf("%d/%d (%.1f%%) TD2 downstream; %d/%d (%.1f%%) TD2 == ref AUG; median downstream offset +%d nt",
                       n_down, nrow(panelC), 100 * n_down / nrow(panelC),
                       n_eq, nrow(panelC), 100 * n_eq / nrow(panelC),
                       if (n_down > 0) as.integer(round(median(panelC$distance_nt[panelC$distance_nt > 0]))) else 0L),
               sprintf("TD2 median %.0f nt vs ref-AUG median %.0f nt (%.1fx)",
                       median(occ_clean$td2_orf_nt),
                       median(occ_clean$ref_atg_orf_nt),
                       median(occ_clean$td2_orf_nt) / median(occ_clean$ref_atg_orf_nt)),
               sprintf("ref median %.2f vs TD2 median %.2f, paired Wilcoxon p=%.3g",
                       median(ref_pwm_occ[koz_valid_occ]),
                       median(td2_pwm_occ[koz_valid_occ]), wt_occ$p.value),
               sprintf("%d/%d (%.1f%%) TD2 downstream of ref AUG; median downstream offset +%d nt",
                       n_down_o, nrow(panelF), 100 * n_down_o / nrow(panelF),
                       if (n_down_o > 0) as.integer(round(median(panelF$distance_nt[panelF$distance_nt > 0]))) else 0L)))
fwrite(summary_dt, file.path(SUPP_DIR, "data", "summary.tsv"), sep = "\t")
cat(sprintf("\nSaved supplement TSVs under %s/data/\n", SUPP_DIR))
