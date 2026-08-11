# Reproducing the analyses

This repository reproduces the analyses in Leshem et al. starting from the Zenodo source-data
record. Raw reads are not required.

## 1. Get the inputs

| What | Where |
|---|---|
| Source data (15 downloads, 695 MB compressed, 4.65 GB extracted) | Zenodo **`10.5281/zenodo.21544336`** — the concept DOI, which always resolves to the latest version. Follow this one, not a version DOI: a version DOI keeps pointing at the files of the version that minted it |
| Tan et al. supplementary tables (§2 only) | already in the record, `tan_2025_supplementary/` |
| MSigDB hallmark gene sets (§3 only) | already in the record, `annotation/msigdb_hallmark_H.tsv` |

Every input is in the record and nothing is fetched at run time.

Unpack so that `source_data/` sits inside the record, and verify **from inside `source_data/`** —
the manifest's paths are relative to it, not to the record root:

```bash
cd source_data && shasum -a 256 -c ../MANIFEST.sha256
```

All 58 entries should report `OK`. Run from the record root instead and it reports
`58 listed files could not be read`, which means the wrong working directory, not a bad download.

## 2. Build the environment

Two environments, because the work ran in two places: the model on a GPU cluster, the figures and
reports on a laptop.

**R** — from the repository root:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore()'
```

[`renv.lock`](renv.lock) pins 187 packages against **R 4.5.2** and **Bioconductor 3.22**, including
`Isopair`, which is not on CRAN — its record points at GitHub `peter4244/Isopair` at commit
`ea33f74`, so `renv::restore()` fetches it without a separate step.

**Python** — one environment per half:

```bash
conda env create -f environment-model.yml     # the deep-learning model
conda env create -f environment-figures.yml   # figures, reports, pandoc
```

Four things these will not do for you:

- **Compilers.** Several Bioconductor packages build from source on Linux. Without a toolchain
  (`build-essential` or equivalent) `renv::restore()` fails partway with a compiler error, not a
  missing-package error. macOS needs the Xcode command-line tools.
- **pandoc.** It is not an R package, so it is not in `renv.lock`; it is in
  `environment-figures.yml`. The reports are rmarkdown and will not render without it.
- **A GPU.** `environment-model.yml` installs `torch==2.5.1+cu121`. With no NVIDIA GPU, replace the
  two `pip:` lines with plain `torch==2.5.1`. CUDA changes speed and floating-point association,
  not the model.
- **Byte-identical figures.** The committed PNGs were rendered under matplotlib **3.10.8**, which
  `environment-figures.yml` pins. See *Known gaps* — fonts are a second and larger reason images
  will differ.

Alternatively, use the container `nmd_1.2.sif` from the same Zenodo record, or build it from the
`Dockerfile` here. See [`ENVIRONMENT.md`](ENVIRONMENT.md).

## 3. Point the repository at the data

Every path resolves through [`config/paths.yml`](config/paths.yml); nothing is hard-coded.

```bash
ln -s /path/to/nmd_deposit_2026 data_deposit
```

Resolution order per key: environment variable `NMD_<KEY>` → `config/paths.yml` → relative to the
repository root. A one-off override is `NMD_DEPOSIT=/elsewhere/source_data Rscript analysis/…`.

Check the wiring before running anything:

```bash
Rscript -e 'source("R/load_config.R"); str(nmd_paths())'
```

`DEPOSIT`, `SQANTI`, `ISOCALL` and `ANNOT` must resolve into the deposit. `CACHE` and `OUT` are
outputs, created on demand. The `LEGACY`/`FIGDATA`/`MASHR_*`/`ISOPAIR` roots point at intermediates
on the authoring machine; they are not needed to reproduce from the deposit and will not resolve
elsewhere.

## 4. Which count matrix?

The deposit ships two long-read isoform matrices because the published analyses used two isoform
universes. They differ by a single isoform, and pairing the wrong matrix with the wrong GTF changes
results silently.

| Section | Counts | Annotation |
|---|---|---|
| §1–2 gene- and isoform-level DE | `nmd_lungcells_counts_4ct.csv` | `sqanti/nmd_lungcells_filtered.gtf` |
| §3–4 isoform-pair / PTC | `nmd_isocall_counts_4ct.csv` | `nmd_isocall_4ct.gtf.gz` |

The shipped code already selects correctly.

## 5. Run

Run the steps in the order below. Outputs go to `tmp/out/` and figures to `figures/`; nothing
writes into the input tree. Wall times are from one laptop run and are there so a twenty-minute
step is not mistaken for a hang.

An `.Rmd` step must be **rendered, not sourced** — `Rscript file.Rmd` will not work.

### 5.1 §1–§2, the report layer (9 steps, ~24 min)

```bash
# NOT named `R` — that would shadow the R binary for the rest of the shell session.
render_rmd() { Rscript -e "rmarkdown::render('$1', quiet=FALSE)"; }

render_rmd analysis/upstream/Isoform_Level_Quantification.Rmd          # 307 s
render_rmd analysis/upstream/NMD_shortread_dge_fullmodel_2026.5.5.Rmd  #  39 s
Rscript    analysis/isopair/isopair_wrapper/01_prepare_data_mashr.R    #  46 s — expression universe
render_rmd analysis/upstream/Isoform_Landscape.Rmd                     #  41 s
render_rmd analysis/upstream/comparison_analysis.Rmd                   #   8 s
render_rmd analysis/upstream/correlation_analysis.Rmd                  #  15 s
render_rmd analysis/upstream/ct_de.Rmd                                 # 314 s
render_rmd analysis/upstream/transcriptional_output.Rmd                # 314 s
render_rmd analysis/upstream/rbp_sr.Rmd                                # 350 s
```

`01_prepare_data_mashr.R` belongs at that position, interleaved rather than grouped with the
isopair steps. Cell-type scope is asserted rather than filtered: the deposited matrices contain
exactly the 26 manuscript samples, and the script stops if the input carries anything else.

### 5.2 §3–§4, the isopair chain (4 steps, ~45 min)

Order matters — `00b` runs *between* `01b` and `02`.

```bash
Rscript analysis/isopair/isopair_wrapper/01b_build_isoform_infrastructure.R  #  949 s — structures, CDS, PTC
Rscript analysis/isopair/isopair_wrapper/00b_build_cds_exons.R              #   24 s — CDS exon table
Rscript analysis/isopair/isopair_wrapper/02_build_profiles_mashr.R          # 1350 s — profiles + pairs
Rscript analysis/isopair/isopair_wrapper/03b_rebuild_cache.R                #  378 s — analysis cache
```

### 5.3 The feature layer (8 steps, ~31 min)

The last step, `export_rds.R`, lives in the model repository (see §5.7) and writes the eight feature
tables the model reads.

```bash
Rscript analysis/isopair/isopair_wrapper/05k_utr5_all_isoforms.R   # 102 s
Rscript analysis/isopair/isopair_wrapper/05s_orfik_scan.R          # 961 s — the long one
Rscript analysis/isopair/isopair_wrapper/05t_ref_cds_features.R    # 209 s
Rscript analysis/isopair/isopair_wrapper/05u_paralog_annotation.R  #   5 s — queries Ensembl live on a cold cache
Rscript analysis/isopair/isopair_wrapper/05s_b_orfik_scan_extend.R # 490 s
Rscript analysis/isopair/isopair_wrapper/05r_ref_atg_analysis.R    #  49 s
Rscript analysis/isopair/isopair_wrapper/05k_b_utr5_refaug.R       #  37 s

# All three settings are required. export_rds.R resolves orfik_scan.rds from $NMD_ISOPAIR_CACHE,
# not from its own directory, and its --results-dir defaults to results_4ct, the PUBLISHED tree.
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R --results-dir results_deposit_h5_2026-08-04   # 9 s
```

### 5.4 The figures the §4 report embeds (~4 min)

The §4 report calls `include_graphics()` on 14 PNGs and **halts** if one is missing, so build them
first. Each panel script reads the exports written by `data_export.R` in its own directory.

```bash
# per-figure data exports
Rscript figures/multipanel/figure3_isopair_and_ptc/data_export.R
Rscript figures/multipanel/figure3_isopair_and_ptc/panel_e_compute.R
Rscript figures/multipanel/figure4_ptcneg_and_model/data_export.R
Rscript figures/supplemental/SF33_CdsAnd3UTR_GENCODE/data_export.R
Rscript figures/supplemental/SF34_TD2Bias_broad/data_export.R

# the panels the report embeds
for f in figures/multipanel/figure3_isopair_and_ptc/figure3_panel*.py \
         figures/multipanel/figure4_ptcneg_and_model/figure4_panel*.py; do python3 "$f"; done

# the composites — Figure 3 and Figure 4 themselves
python3 figures/multipanel/figure3_isopair_and_ptc/figure3_composite.py
python3 figures/multipanel/figure4_ptcneg_and_model/figure4_composite.py
```

### 5.5 The §3 and §4 reports

```bash
render_rmd analysis/upstream/productive_response.Rmd                                      # 855 s — §3
render_rmd analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd  #  42 s — §4
```

### 5.6 Supplemental figures for §1–§4

Fourteen producers, run by one generator. **Run it from that directory** — the scripts find their
inputs through `config/paths.yml`, but write output relative to the working directory.

```bash
cd analysis/upstream/Figures
Rscript -e 'rmarkdown::render("make_supplemental_figures.Rmd")'
```

It needs `structures.rds` and `cds.rds` from step 5.2, and the DGELists and per-cell-type mashr
outputs under `tmp/out/mashr_gene/` and `tmp/out/mashr_isoform/` from step 5.1, so run those first.
`_nmd_machinery_report.R` sits in the same directory and is deliberately not run.

These rebuild the supplement from the same data as the main text; they are not a re-check of the
published supplemental figures, and the two differ.

### 5.7 §5, the deep-learning model

§5 is the one section that does not run from this repository alone. It needs the model repository
**`NMD_orf_model_v5_4ct`** ([10.5281/zenodo.21536501](https://doi.org/10.5281/zenodo.21536501)) and
a GPU. Clone it beside this one; `<model repo>` below is its path, and the `python` steps run from
inside it.

**You do not need to run any of this to reproduce §5's reported numbers.** The deposit's
`model.zip` carries the trained checkpoint, the per-isoform predictions and the full
interpretability export, and the §5 figure scripts read those. This chain rebuilds them.

**Submit the wrappers rather than typing the commands.** The model repository carries a
`slurm_*_dn.sh` driver per stage (`_dn` = deposit-native). They default `RESULTS_DIR` to the
deposit-native tree and derive the window from `paths_config.py --selected-tag` instead of
repeating a literal. The command each runs is shown so the page can be read without a cluster.

**Three result directories, two of them traps:**

| directory | what it is |
|---|---|
| `results_4ct` | the **published** run, and the default of every script taking `--results-dir`. Writing here destroys what the rebuild is compared against |
| `results_4ct_dn` | **deprecated** — its HDF5 was built from Channing inputs; now under `deprecated_2026-08-04/` |
| `results_deposit_h5_2026-08-04` | **the deposit-native tree**, and what every `_dn` wrapper defaults to |

**The window is 1000/1000 and the config default is not.** The deposited checkpoint is
`best_model_atg1000_stop1000_seed42.pt`. `config_dn.yaml` carries the window twice —
`window_size_atg: 1000` under `selected:` and `window_size_atg: 100` under `data:`, the sweep grid's
starting point — and `03_train.py` reads the second. A run that does not pass the window trains at
100/1000. Anything you write yourself must pass `--atg-window 1000 --stop-window 1000`.
`window_size` is a full width: `data_prep.py` takes `half_win = win_size // 2`, so 1000 means ±500.

The model does not take bare sequence: `data_prep.py`'s nine channels include exon-junction
positions and the structural features include `n_downstream_ejc`, so scoring an arbitrary RNA
sequence is not a well-posed use of it.

**1 — the feature tables** (local, no GPU). Already run as the last step of 5.3; repeated here
because the rest of this chain depends on it.

```bash
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R --results-dir results_deposit_h5_2026-08-04
```

**2 — build the HDF5.** `sbatch slurm_build_h5_dn.sh`

```bash
python data_prep.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 --workers 8
```

**3 — patch the stop-codon column.** `sbatch slurm_patch_selected_orfs.sh` — but note that
`patch_stop_codon.py` defaults `--results-dir` to the published tree and the wrapper passes no
arguments, and that without `--in-place` it writes `selected_orfs_fixed.tsv` beside the original.
Required for §4.1's χ² test.

```bash
python scripts/patch_stop_codon.py --results-dir results_deposit_h5_2026-08-04 --in-place
```

**4 — the determinism gate.** Must pass before training, at the real window config (~30 s).

```bash
python verify_determinism.py --atg 1000 --stop 1000
```

**5 — train.** `sbatch slurm_train_dn.sh`, which also scores the `val_clean` development split.

```bash
python 03_train.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000
python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000 --member-seed 42 --split val_clean
```

**6 — the interpretation layer**, which runs over the full cohort, not the test split.

```bash
sbatch slurm_uorf_dn.sh              # infer_uorf_attention.py, then the uORF metrics
sbatch slurm_deepshap_joint_dn.sh    # deepshap.py --branches joint, 5 runs
sbatch slurm_deepshap_structural_dn.sh
sbatch slurm_kernel_shap_dn.sh       # 11_kernel_shap_branches.py
sbatch slurm_export_chain_dn.sh      # 06_, 07_, 08_, 09_ GC/junction/polyA, 09b_, 09c_, 09d_
```

`slurm_deepshap_all_dn.sh` does joint, structural and `atg stop` in one job. The five DeepSHAP runs
vary the background seed; the replicate mean and sd are the uncertainty statement.

**7 — the final evaluation. Run once, and last.**

```bash
sbatch slurm_eval_final_dn.sh
# python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
#     --atg-window 1000 --stop-window 1000 --member-seed 42 --split test_clean --final
```

The test split is untouched by everything above: training early-stops on `val_clean` and the
interpretation chain runs `--split all --full-cohort`. `--final` marks the single evaluation
allowed to touch chr1/3/5/7, and the metrics JSON records `evaluation_class=final_test`.

`test_clean` is the paralog-free held-out set — chr1/3/5/7, with the 122 paralog-straddling
transcripts under their own `test_paralog` label. It is the split the manuscript describes;
`--split test` is a different population and gives different numbers.

**8 — the report and the §5 figure data.**

```bash
NMD_RESULTS_DIR=results_deposit_h5_2026-08-04 \
  Rscript -e 'rmarkdown::render("<model repo>/orf_model_report_v5.Rmd",
                                knit_root_dir = normalizePath("<model repo>"))'

# Both variables are required. MODEL_RESULTS is the model's own results; FEATURES is what
# export_rds.R wrote with --results-dir. SF40's data_export.R needs the same one.
NMD_MODEL_RESULTS=<model repo>/results_deposit_h5_2026-08-04 \
NMD_FEATURES=<model repo>/results_deposit_h5_2026-08-04 \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R
```

**Flags whose omission fails silently or late:**

- **`--member-seed` is required by every consumer of the checkpoint.** `03_train.py` writes
  `best_model_{tag}_seed{N}.pt`, one file per ensemble member, so `evaluate.py`,
  `11_kernel_shap_branches.py` and `deepshap.py` refuse to load an unqualified
  `best_model_{tag}.pt`. Omit it and the step dies with `FileNotFoundError` *after* training has
  succeeded. `--member-seed` names the checkpoint; `--seed` is the RNG seed.
- **`--results-dir` is not optional, and not only for the training scripts.** `export_rds.R`,
  `infer_uorf_attention.py`, `compute_uorf_attention_metrics.R` and `orf_model_report_v5.Rmd` all
  default to the published run, so omitting it reads or writes the wrong vintage silently.
- **The report takes `NMD_RESULTS_DIR`, not a flag** — `rmarkdown::render()` does not pass
  arguments through to the document.
- **DeepSHAP cannot run deterministically.** `shap`'s DeepLIFT routes MaxPool gradients through
  `max_unpool1d`, which has no deterministic kernel, so it needs `NMD_ALLOW_NONDETERMINISM=1`.
  The five replicates are the uncertainty statement instead.
- **Joint and structural DeepSHAP write the same filename** and will overwrite each other, or be
  globbed and averaged together by the consumers. `--branches joint` is what Figure 5D used.
- **`relabel_tx_summary_4ct.R` is not needed** — a proven no-op on a deposit-native scaffold.

Determinism is enforced but bounded: `set_seed` pins `cudnn.deterministic`, `CUBLAS_WORKSPACE_CONFIG`
and `torch.use_deterministic_algorithms`, and two seeded runs are then bit-identical *on the same
GPU model and library versions* (verified on a Tesla V100-SXM2-32GB, torch 2.5.1+cu121). A pass does
not cover the GPU training path, which runs under AMP, and `NMD_ALLOW_NONDETERMINISM=1` disables all
of it without leaving a trace in the outputs.

What the deposited run produces, read from the deposit's own
`model/metrics_atg1000_stop1000_seed42_test_clean.json`: **AUC 0.9257, AUPRC 0.8175, n_eval 10,522
with 2,405 NMD susceptible**, `evaluation_class = final_test`, `best_epoch = 5`, window 1000/1000,
`member_seed = 42`. These are the numbers §5 reports at its own precision, 0.93 and 0.82.

## What is verified, and what is not

[`docs/VERIFIED.md`](docs/VERIFIED.md) records, per claim, what the rebuild has been observed to do,
including which claims carry a value that differs from the published one. It is generated from the
claims ledger, not written by hand.

## Known gaps

- **Set `NMD_CLAIM_VALUES` before you start, or lose the claim table.** Ten of the `.Rmd` steps call
  `claim_emit()`, which writes one row per reported quantity — value, n, population, producer — to
  `$NMD_CLAIM_VALUES`, or to `tmp/claim_values.tsv` when unset. `tmp/` is git-ignored, so the file
  is written and then lost with the working directory. It is the machine-readable record of what
  your run computed:

  ```bash
  export NMD_CLAIM_VALUES="$PWD/logs/claim_values.tsv"
  ```

- **Your figures will not be byte-identical to the paper's.** The published figures were rendered on
  macOS in **Arial**. The container has no Arial — it is licensed and cannot be redistributed — so
  it renders **Liberation Sans**, a metric-compatible clone. Advance widths and therefore layout are
  identical; glyph outlines are not, so a pixel diff differs everywhere while panel positions and
  wrapping match. Figures are therefore verified by comparing their **data exports**, not their
  images. The font order is set once in `figures/lib/ggplot_style.py`:
  `Arial → Helvetica Neue → Helvetica → Liberation Sans → Nimbus Sans → DejaVu Sans`.

- **Stages beyond Isopair stage 1 are wired but not verified end-to-end.** Which claims that covers
  is in `docs/VERIFIED.md` under *Not measured*; an unmeasured claim is not an unaffected one.

- **The §5 chain needs a GPU and cluster access**, so it is the one part of this document that
  cannot be run from the deposit alone on a laptop.
