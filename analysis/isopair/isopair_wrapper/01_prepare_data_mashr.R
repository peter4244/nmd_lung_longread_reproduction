#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Data Preparation (mashr Classification)
# Version: 7.0-mashr
#
# Purpose: Load isocall count matrix, normalize, filter, classify isoforms as
#          NMD-sensitive or non-NMD using mashr differential isoform expression
#          results. All study-specific logic lives here.
#
# Input:
#   - Isocall count matrix + GTF
#   - 6 mashr DIE CSVs (one per cell type)
#
# Output (to output_dir/):
#   - expression_data.rds    — filtered CPM matrix (isoform x sample)
#   - sample_metadata.rds    — data.frame with sample_id, treatment, ct, donor
#   - gene_map.rds           — tibble with isoform_id, gene_id
#   - nmd_classification.rds — named list: per cell type + all_samples
#   - dmso_samples.rds       — named list of DMSO sample IDs per cell type
#   - smg1i_samples.rds      — named list of Smg1i sample IDs per cell type
#
# Usage:
#   Rscript 01_prepare_data_mashr.R [--output-dir DIR]

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
library(edgeR)

# ==============================================================================
# 0. Configuration
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
# Config-resolved so a stray cwd cannot overwrite the retained reference intermediates
# in the legacy tree. Override with --output-dir.
output_dir <- file.path(.P$OUT, "data_mashr")
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  if (idx < length(args)) output_dir <- args[idx + 1]
}
# Rule in code: never write into a tree we read as input (R/load_config.R).
# MISSED in the 2026-07-26 sweep -- this script takes --output-dir and writes six
# reference objects, and was the only writer left unguarded. Found by review.
nmd_assert_writable(output_dir, .P, "--output-dir")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Input paths ---
count_matrix_file <- file.path(.P$DEPOSIT, "nmd_isocall_counts_4ct.csv")
gtf_file <- file.path(.P$DEPOSIT, "nmd_isocall_4ct.gtf.gz")

# mashr DE results.
#
# REGENERATION-FIRST (2026-07-26, Pete's decision: documented regeneration rather than
# depositing the four 57 MB CSVs). The producer is
# Isoform_Level_Quantification.Rmd :: mashr-save, which reads deposited inputs only and
# derives its output suffix from the mashr COEFFICIENT names -- Smg1i_in_AT2, Smg1i_in_LAE,
# ... -- so a regenerated file is nmd_mashr_die_at2_*, _lae_*. The retained copies on disk
# predate the AT->AT2 / DD->LAE rename and are _at_, _dd_.
#
# Asking for one naming and finding the other fails loudly at read.csv, which is the good
# case; the bad case is a consumer that globs and silently picks up either vintage. So this
# resolves EXPLICITLY, prefers the regenerated file, and PRINTS which vintage it used --
# the vintage is provenance and must appear in the log, not be inferred later.
de_date <- "2026.3.10"
# Paper scope (2026-04-29): 4 non-ALI cell types only — AT, DD, FB, MV (internal codes).
# DD_ALI / DO_ALI / DO dropped from the paper.
de_suffix_current <- c(AT = "at2", DD = "lae", FB = "fb", MV = "mv")   # post-rename
de_suffix_legacy  <- c(AT = "at",  DD = "dd",  FB = "fb", MV = "mv")   # pre-rename, retained
# Regenerated output goes to .P$OUT/mashr_isoform -- that is where the producer actually
# writes (Isoform_Level_Quantification.Rmd:591, OUT_DIR_MASH <- file.path(.P$OUT,
# "mashr_isoform")). The first version of this searched .P$OUT/mashr, which the producer
# never writes, and that one wrong word made the whole guard useless: see below.
# THE DEPOSIT IS SEARCHED FIRST, added 2026-08-10. These four DIE tables were reachable only
# from a regeneration directory or, behind NMD_ALLOW_LEGACY_MASHR=1, from the pre-migration
# tree -- neither of which an outside reproducer has. They are now in the Zenodo record
# (source_data/isopair/, 57 MB), so a reader with the deposit can run this script.
de_regen_dirs <- c(file.path(.P$DEPOSIT, "isopair"),
                   file.path(.P$OUT, "mashr_isoform"), file.path(.P$OUT, "mashr"))

# THE LEGACY FALLBACK IS OPT-IN AND DEFAULTS OFF, and that is the important part.
#
# A fallback that always succeeds converts a loud failure into a silent one. With the
# retained CSVs sitting in .P$MASHR_LR, the stop() below could never fire: a regeneration
# writing to an unsearched directory would be skipped, the legacy vintage returned, and a
# filename printed that reads as success. The re-measurement would then hand back the
# PUBLISHED 1,548 pairs / 51.0% SE and look like proof that regeneration reproduces them.
# Vacuous, and confidently so. For a proof run this must be off, so absence fails loudly.
allow_legacy <- identical(Sys.getenv("NMD_ALLOW_LEGACY_MASHR"), "1")
de_search_dirs <- if (allow_legacy) c(de_regen_dirs, .P$MASHR_LR) else de_regen_dirs

resolve_de_file <- function(ct) {
  for (d in de_search_dirs) for (sfx in list(de_suffix_current, de_suffix_legacy)) {
    f <- file.path(d, sprintf("nmd_mashr_die_%s_%s.csv", sfx[[ct]], de_date))
    if (file.exists(f)) return(f)
  }
  stop(sprintf(paste0("no mashr DIE file for %s. Looked for suffixes {%s, %s} under {%s}.\n",
                      "  Regenerate via Isoform_Level_Quantification.Rmd :: mashr-save, which\n",
                      "  writes to .P$OUT/mashr_isoform. To fall back to the retained\n",
                      "  pre-rename copies instead, set NMD_ALLOW_LEGACY_MASHR=1 -- but NOT\n",
                      "  during a regeneration proof run, or the result is vacuous."),
               ct, de_suffix_current[[ct]], de_suffix_legacy[[ct]],
               paste(de_search_dirs, collapse = ", ")))
}

# --- Thresholds ---
# mashr has pre-computed nmd_responsive column; non-NMD uses adj.P.Val.
# Threshold raised from 0.50 → 0.30 (2026-04-29) because the 4-CT mashr rerun
# tightened adj.P.Val distributions (the per-CT limma fits feeding mashr were
# rerun on the 4-CT subset). At 0.50 only ~7-12% of isoforms qualified as
# non-NMD vs ~71% historically; 0.30 yields ~52-60% non-NMD per CT, restoring
# a usable pair set for downstream gene-matching.
NON_NMD_ADJ_P <- 0.30
MIN_PROP      <- 0.05

# PLATFORM-DETERMINISM TOLERANCE ON THE 5% COMPARISON (2026-08-10).
#
# `max_prop >= MIN_PROP` is not platform-deterministic, and the consequence reaches the paper.
# Measured on ONE machine, same script, same deposited counts, same edgeR 4.8.2: macOS R passes
# 207,913 isoforms here and the Linux container passes 207,987, ending at 95,623 and 95,639
# expressed isoforms and therefore at 1,578 vs 1,580 gene-matched pairs.
#
# The reason is that max_prop is a per-gene share, so a gene with twenty comparable isoforms puts
# its members ON 0.05 rather than near it. 1,234 isoforms cross this comparison between the two
# platforms; 878 of them compute to EXACTLY 0.05 on one side, and the largest distance to the
# threshold across all 1,234 is 3.5e-17. The largest disagreement in max_prop anywhere is 1.8e-15.
# The platforms are not computing different quantities -- they are ordering the same sums
# differently, and a value mathematically equal to 0.05 lands a few ulp either side.
#
# So the comparison is taken at a tolerance far above the noise (1e-9) and far below any real
# effect: no isoform's biology places it within a billionth of a 5% share. A share that IS 5%
# passes a ">= 5%" filter, which is also what the filter means.
PROP_TOL      <- 1e-9

cat("=== NMD Data Preparation — mashr Classification ===\n\n")

# ==============================================================================
# 1. Load count matrix and build DGEList
# ==============================================================================

cat("Loading count matrix...\n")
count_df <- read.csv(count_matrix_file, check.names = FALSE)
isoform_ids <- count_df$id
count_mat <- as.matrix(count_df[, -1])
rownames(count_mat) <- isoform_ids
rm(count_df); gc(verbose = FALSE)
cat(sprintf("  %d isoforms x %d samples\n", nrow(count_mat), ncol(count_mat)))

# Parse sample metadata from column names
# Deposit format: <donor>_<treatment>_<cellType>, e.g. 001V_DMSO_AT2.
# The published cell-type names (AT2/LAE) are mapped back to this pipeline's internal
# codes (AT/DD) so every downstream object keeps the key space it has always used.
CT_FROM_PUBLISHED <- c(AT2 = "AT", LAE = "DD", FB = "FB", MV = "MV")
parse_sample <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  stopifnot("unexpected deposit column name" = length(parts) == 3)
  ct <- CT_FROM_PUBLISHED[[parts[3]]]
  data.frame(sample_id = x, ct = ct, donor = parts[1], treatment = parts[2],
             stringsAsFactors = FALSE)
}

sample_metadata <- do.call(rbind, lapply(colnames(count_mat), parse_sample))
cat(sprintf("  %d samples: %s\n", nrow(sample_metadata),
            paste(names(table(sample_metadata$ct)), collapse = ", ")))

# --- 4-CT-native input: nothing is dropped here ---
# The deposited matrices contain exactly the 26 manuscript samples (13 donor-matched
# DMSO/SMG1i pairs across AT/DD/FB/MV). Two exclusions that earlier versions of this
# script performed are therefore already applied at source, and are asserted rather
# than re-applied:
#
#   1. DO donor 029T, a PCA outlier (DMSO sample dist=252 from the DO centroid, >3x
#      the next), both treatment partners removed to keep the design balanced.
#   2. DD_ALI and DO, out of manuscript scope.
#
# Asserting instead of filtering is the point of the 4-CT-native rewrite: the isoform
# universe, TMM normalization, the 5% and filterByExpr filters, reference selection,
# the 25% floor and the C2 NMD partner are AT/DD/FB/MV by construction, with zero
# 6-CT trace. If the input ever regains out-of-scope samples this fails loudly here
# rather than silently changing every downstream result.
# See REFERENCE_FLOOR_PLAN.md AMENDMENT 2026-07-11 (decision b), R1/M-1.
manuscript_cts <- c("AT", "DD", "FB", "MV")
stopifnot(
  "input carries cell types outside the manuscript scope" =
    setequal(unique(sample_metadata$ct), manuscript_cts),
  "input is not the expected 26 manuscript samples" = ncol(count_mat) == 26L,
  "count_mat and sample_metadata are misaligned" =
    identical(colnames(count_mat), sample_metadata$sample_id),
  "expected 13 donor-matched DMSO/SMG1i pairs" =
    all(table(sample_metadata$treatment) == 13L)
)
cat(sprintf("  4-CT-native input verified (%s): %d samples, %d donors\n",
            paste(manuscript_cts, collapse = ", "), ncol(count_mat),
            length(unique(sample_metadata$donor))))

# Build gene map from GTF
cat("Parsing gene map from GTF...\n")
gtf_gr <- rtracklayer::import(gtf_file)
tx_rows <- as.data.frame(gtf_gr[gtf_gr$type == "transcript", ])
gene_map <- data.frame(
  isoform_id = tx_rows$transcript_id,
  gene_id    = tx_rows$gene_id,
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!duplicated(gene_map$isoform_id), ]
cat(sprintf("  %d isoforms mapped to %d genes\n",
            nrow(gene_map), length(unique(gene_map$gene_id))))

# DGEList + TMM normalization
dge <- DGEList(counts = count_mat,
               genes = data.frame(gene_id = gene_map$gene_id[
                 match(rownames(count_mat), gene_map$isoform_id)]))
dge <- calcNormFactors(dge, method = "TMM")
cpm_mat <- cpm(dge, log = FALSE)

# ==============================================================================
# 2. Condition-stratified 5% isoform filter
# ==============================================================================

cat("\nApplying 5% condition-stratified filter...\n")
dmso_cols  <- sample_metadata$sample_id[sample_metadata$treatment == "DMSO"]
smg1i_cols <- sample_metadata$sample_id[sample_metadata$treatment == "Smg1i"]

# Per-gene proportions in each condition
compute_max_prop <- function(cpm, samples, gene_ids) {
  cpm_sub <- cpm[, samples, drop = FALSE]
  # Gene totals per sample
  gene_totals <- rowsum(cpm_sub, gene_ids)
  # Isoform proportions: cpm / gene_total for its gene
  props <- cpm_sub / gene_totals[gene_ids, , drop = FALSE]
  props[is.nan(props)] <- 0
  # Max proportion across samples within this condition
  apply(props, 1, max, na.rm = TRUE)
}

gene_ids_all <- gene_map$gene_id[match(rownames(cpm_mat), gene_map$isoform_id)]
max_prop_dmso  <- compute_max_prop(cpm_mat, dmso_cols, gene_ids_all)
max_prop_smg1i <- compute_max_prop(cpm_mat, smg1i_cols, gene_ids_all)

keep_5pct <- max_prop_dmso >= (MIN_PROP - PROP_TOL) | max_prop_smg1i >= (MIN_PROP - PROP_TOL)
cat(sprintf("  5%% filter: %d / %d isoforms pass\n", sum(keep_5pct), length(keep_5pct)))

# filterByExpr with design (cell type + treatment)
design <- model.matrix(~ 0 + ct + treatment, data = sample_metadata)
keep_expr <- filterByExpr(dge[keep_5pct, ], design = design,
                          min.count = 2, min.total.count = 4)
keep_final <- names(keep_5pct)[keep_5pct]
keep_final <- keep_final[keep_expr]

# ---- THE DEPOSITED UNIVERSE WINS, AND THE RECOMPUTATION IS THE FALLBACK ---------------------
#
# The filter above is now platform-deterministic (see PROP_TOL) but it does NOT reproduce the
# universe the manuscript was written from, and no threshold rule can. 2,567 isoforms lie within
# 1e-9 of the 5% boundary, and the published 95,623 is one particular floating-point split of that
# band -- the split macOS happened to make in July. Recomputing here, deterministically, gives
# 95,817 and 1,580 / 5,052 pairs against the paper's 1,578 / 5,049.
#
# So the universe is DEPOSITED as data rather than re-derived as a result. A reader who reproduces
# from the Zenodo record reads the same 95,623 isoforms this paper used and gets its numbers
# exactly, on any platform. A reader regenerating from counts alone gets the deterministic
# recomputation and a universe that differs by ~0.2%; REPRODUCTION.md states both and why.
#
# WHICH PATH RAN IS PRINTED, ALWAYS. A pipeline that silently prefers a deposited answer to a
# computed one is indistinguishable from a pipeline that computed it, and that is the whole failure
# this repository keeps paying for.
universe_file <- file.path(.P$DEPOSIT, "isopair", "expressed_isoforms_4ct.txt")
if (file.exists(universe_file)) {
  deposited <- readLines(universe_file)
  extra   <- setdiff(keep_final, deposited)
  missing <- setdiff(deposited, keep_final)
  cat(sprintf("  UNIVERSE: deposited (%s)\n", basename(universe_file)))
  cat(sprintf("            %d isoforms; recomputation here gave %d (%d not in the deposit, %d only in it)\n",
              length(deposited), length(keep_final), length(extra), length(missing)))
  # The deposited universe must be a set this pipeline could have produced. Anything in it that the
  # count matrix does not contain means the two are from different vintages, which is a stop.
  absent <- setdiff(deposited, rownames(cpm_mat))
  if (length(absent))
    stop(sprintf("deposited universe names %d isoform(s) absent from the count matrix, e.g. %s -- the deposit and the counts are different vintages",
                 length(absent), paste(utils::head(absent, 3), collapse = ", ")))
  keep_final <- deposited
} else {
  cat(sprintf("  UNIVERSE: RECOMPUTED -- %s not found.\n", universe_file))
  cat("            This does not reproduce the manuscript's 95,623; see REPRODUCTION.md.\n")
}

# ---- OPEN ISSUE (W42): FUSION GENES ARE RETAINED, DELIBERATELY --------------------
#
# There used to be an exclusion here:
#
#     fusion_genes    <- unique(gene_map$gene_id[grepl("--", gene_map$gene_id)])
#     fusion_isoforms <- gene_map$isoform_id[gene_map$gene_id %in% fusion_genes]
#     keep_final      <- setdiff(keep_final, fusion_isoforms)
#
# IT NEVER REMOVED ANYTHING. gene_map$gene_id comes from the GTF, and the deposited
# nmd_isocall_4ct.gtf.gz separates fusion partners with "::" -- 671 distinct fusion gene_ids
# in a 400k-transcript sample, ZERO containing "--". So the grep matched nothing, the setdiff
# was a no-op, and the cat() below announced a filter that had not run. A silent zero.
#
# It was in effect for the PUBLISHED run too: the published universe is 162,800 and the
# rebuild reproduces 162,800 exactly, which it could not if the exclusion had ever fired.
#
# RETAINED ON PURPOSE (Pete, 2026-07-26), because the manuscript already documents them as
# present. §1's structural-category sentence reports "1.9% falling into other transcript
# classes, including fusion and genic isoforms", and that 1.9% is fusion 2,719 + genic 406 =
# 3,125 -- verified EXACT (claim 1.1.3). The paper and the data agree; only the old comment
# disagreed.
#
# Measured cost of excluding them instead (gene_id containing "::", the intended definition):
#   162,800 -> 162,365 expressed isoforms   |  34,387 -> 34,269 NMD susceptible
#   7,105   -> 7,073   NMD susceptible genes|  per cell type -66 to -110
# i.e. every affected headline moves by under 0.5% and no conclusion changes -- which is why
# retaining is preferred over restating four abstract numbers.
#
# STILL OPEN, and the reason this stays an issue rather than a closed decision:
#   (a) "fusion" has two incompatible definitions here -- gene_id "::" selects 435 isoforms,
#       SQANTI structural_category == "fusion" selects 2,719, and only 330 satisfy both. The
#       retained set is whichever the reader assumes, and the Methods does not say.
#   (b) fusion gene_ids being present is what makes the greedy sub("\\..*","") collapse at
#       five aggregation sites live rather than hypothetical (W41).
# Do not re-add an exclusion without settling (a).

cpm_mat <- cpm_mat[keep_final, ]
cat(sprintf("  After filterByExpr: %d isoforms (fusion genes RETAINED -- see W42)\n",
            nrow(cpm_mat)))

# ==============================================================================
# 3. NMD classification per cell type (mashr)
# ==============================================================================

cat("\nClassifying isoforms using mashr results...\n")
nmd_classification <- list()
dmso_samples_list  <- list()
smg1i_samples_list <- list()

for (ct in names(de_suffix_current)) {
  de_file <- resolve_de_file(ct)
  cat(sprintf("  %s: mashr DIE <- %s\n", ct, de_file))
  de <- read.csv(de_file)

  # mashr uses 'txid' column (not 'transcript_id')
  # Filter to isoforms in our expression data
  de <- de[de$txid %in% rownames(cpm_mat), ]

  # mashr provides pre-computed nmd_responsive column
  nmd_ids     <- de$txid[de$nmd_responsive == TRUE]
  non_nmd_ids <- de$txid[de$adj.P.Val > NON_NMD_ADJ_P]

  nmd_classification[[ct]] <- list(nmd = nmd_ids, non_nmd = non_nmd_ids)

  ct_samples <- sample_metadata$sample_id[sample_metadata$ct == ct]
  dmso_samples_list[[ct]]  <- ct_samples[ct_samples %in% dmso_cols]
  smg1i_samples_list[[ct]] <- ct_samples[ct_samples %in% smg1i_cols]

  cat(sprintf("  %s: %d NMD, %d non-NMD\n", ct, length(nmd_ids), length(non_nmd_ids)))
}

# all_samples: NMD = union across AT, DD, FB, MV (paper scope: 4 non-ALI cell types)
# non-NMD = intersection across AT, DD, FB, MV
main_cts <- c("AT", "DD", "FB", "MV")
all_nmd <- unique(unlist(lapply(nmd_classification[main_cts], `[[`, "nmd")))
all_non_nmd <- Reduce(intersect, lapply(nmd_classification[main_cts], `[[`, "non_nmd"))
nmd_classification[["all_samples"]] <- list(nmd = all_nmd, non_nmd = all_non_nmd)
dmso_samples_list[["all_samples"]]  <- dmso_cols
smg1i_samples_list[["all_samples"]] <- smg1i_cols
cat(sprintf("  all_samples (AT, DD, FB, MV): %d NMD (union), %d non-NMD (intersection)\n",
            length(all_nmd), length(all_non_nmd)))

# ==============================================================================
# 4. Save outputs
# ==============================================================================

cat(sprintf("\nSaving to %s/...\n", output_dir))
saveRDS(cpm_mat,              file.path(output_dir, "expression_data.rds"))
saveRDS(sample_metadata,      file.path(output_dir, "sample_metadata.rds"))
saveRDS(gene_map,             file.path(output_dir, "gene_map.rds"))
saveRDS(nmd_classification,   file.path(output_dir, "nmd_classification.rds"))
saveRDS(dmso_samples_list,    file.path(output_dir, "dmso_samples.rds"))
saveRDS(smg1i_samples_list,   file.path(output_dir, "smg1i_samples.rds"))

cat("\nDone.\n")
