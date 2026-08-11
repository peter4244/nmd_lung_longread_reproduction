# Reproducing the analyses

This repository reproduces the analyses in Leshem et al. starting from the Zenodo source-data
record. Raw reads are not required.

## 1. Get the inputs

| What | Where |
|---|---|
| Source data (15 downloads, 695 MB compressed, 4.65 GB extracted) | Zenodo **`10.5281/zenodo.21544336`** — the CONCEPT DOI, which always resolves to the latest version. Cite and follow this one, not a version DOI: publishing a new version mints a new record id, and a version DOI silently keeps pointing at the old files |
| Tan et al. supplementary tables (§2 reanalysis only) | **in the Zenodo record** (`tan_2025_supplementary/`) — you do not fetch these separately |
| MSigDB hallmark gene sets (§3 only) | **in the Zenodo record** (`annotation/msigdb_hallmark_H.tsv`) — you do not fetch these separately |

Unpack the record so that `source_data/` sits inside it, and verify **from inside
`source_data/`** — the manifest's paths are relative to it, not to the record root:

```bash
cd source_data && shasum -a 256 -c ../MANIFEST.sha256
```

Every one of the 58 entries should report `OK`. Running it from the record root instead reports
`58 listed files could not be read`, which is a wrong working directory and not a corrupt download.

**Every input is in the Zenodo record; nothing is fetched at run time.** `analysis/upstream/productive_response.Rmd` used to call `msigdbr` at knit time, which downloads MSigDB into a per-user cache and therefore failed on any machine whose cache was cold — including every clean container. The hallmark collection is now deposited as `annotation/msigdb_hallmark_H.tsv` (50 gene sets, 7,322 gene-set/symbol pairs, MSigDB **2026.1.Hs**) and the chunk reads it. It is regenerated, and checked against MSigDB, by `analysis/upstream/capture_msigdb_hallmark.R` — a producer, not a run-time dependency.

## 2. Build the environment

**Two environments, not one, because the work genuinely ran in two places** — the model on a GPU
cluster, the figures and reports on a laptop. Merging them would mean choosing one matplotlib, and
matplotlib is the one package that cannot be chosen freely (below).

**R** — from the repository root:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore()'
```

[`renv.lock`](renv.lock) pins 187 packages against **R 4.5.2** and **Bioconductor 3.22**. It
includes `Isopair`, which is not on CRAN or Bioconductor; its record points at GitHub
`peter4244/Isopair` at commit `ea33f74` (tag `v1.0.0`, archived as `10.5281/zenodo.21536494`), so
`renv::restore()` fetches it without a separate step.

**Python** — one environment per half:

```bash
conda env create -f environment-model.yml     # the deep-learning model
conda env create -f environment-figures.yml   # figures, reports, pandoc
```

Four things this will not do for you, each of which has actually bitten someone here:

- **Compilers.** Several Bioconductor packages have C/C++ sources and build from source on Linux.
  Without a toolchain (`build-essential` or equivalent) `renv::restore()` fails partway with a
  compiler error, not a missing-package error. macOS needs the Xcode command-line tools.
- **pandoc is not an R package.** The reports are knitr/rmarkdown and will not render without it,
  and renv cannot see it, so it does not appear in `renv.lock` — it is carried in
  `environment-figures.yml` instead. On the first run of this pipeline on a second machine, pandoc
  was the very first thing to fail: *naming an interpreter does not bring its toolchain.*
- **A GPU.** `environment-model.yml` installs `torch==2.5.1+cu121` from PyTorch's own index. With
  no NVIDIA GPU, replace the two `pip:` lines with plain `torch==2.5.1`. CUDA changes speed and
  floating-point association, not the model; attribution values are sums of floats and may differ
  in their last decimal places between CPU and GPU.
- **Byte-identical figures.** The committed PNGs were rendered under matplotlib **3.10.8**, which is
  what `environment-figures.yml` pins. A different patch release changes rasterisation, so a PNG
  that does not byte-compare is a rendering difference and not a failure to reproduce — and see the
  font note under *Known gaps*, which is a second and larger reason images will differ.

> **Two honest limits, recorded rather than smoothed over.** `hexbin` and `mclust` are referenced by
> code in this repository but are installed on none of the machines we can inspect, so `renv.lock`
> has no record of them and a restored environment will not contain them. And **these three files have not yet been used to build an
> environment from scratch and run the chain inside it** (`W292`). A lockfile that has never been
> restored is an assertion, and it is labelled as one here until that run happens.

## 3. Point the repository at it

Every path in this repository resolves through [`config/paths.yml`](config/paths.yml). Nothing
is hard-coded.

```bash
ln -s /path/to/nmd_deposit_2026 data_deposit
```

Resolution order for each key is: environment variable `NMD_<KEY>` → the value in
`config/paths.yml` → interpreted relative to the repository root. So a one-off override is just:

```bash
NMD_DEPOSIT=/elsewhere/source_data Rscript analysis/…
```

Check the wiring before running anything:

```bash
Rscript -e 'source("R/load_config.R"); str(nmd_paths())'
```

`DEPOSIT`, `SQANTI`, `ISOCALL` and `ANNOT` must resolve into the deposit. `CACHE` and `OUT` are
outputs and are created on demand. The `LEGACY`/`FIGDATA`/`MASHR_*`/`ISOPAIR` roots point at
retained intermediates on the authoring machine; they are **not** needed to reproduce from the
deposit and will simply not resolve elsewhere.

## 4. Which count matrix?

The deposit ships two long-read isoform matrices because the published analyses used two
different isoform universes. Use the one matching the section you are reproducing — they differ
by a single isoform, and pairing the wrong one with the wrong GTF silently changes results.

| Section | Counts | Annotation |
|---|---|---|
| §1–2 gene- and isoform-level DE | `nmd_lungcells_counts_4ct.csv` | `sqanti/nmd_lungcells_filtered.gtf` |
| §3–4 isoform-pair / PTC | `nmd_isocall_counts_4ct.csv` | `nmd_isocall_4ct.gtf.gz` |

The deposit README explains the difference. The shipped code already selects correctly.

## 5. Run

Outputs are written under `tmp/out/` (`.P$OUT`) and figures under `figures/`. Nothing writes
into the input tree.

**This sequence is the one that EXECUTED, not a reconstruction.** Stages 1–3 ran 9/9, 4/4 and 8/8
at exit 0 in a throwaway projection reading only the deposit (C68, C70, C72), and §3/§4's two
reports ran at exit 0 afterwards (C74). The wall times are from those runs, on one laptop, and are
there so a step that takes twenty minutes is not mistaken for a hang.

**An `.Rmd` step must be rendered, not sourced.** The executed runs wrapped each one in
`rmarkdown::render()`; `Rscript file.Rmd` will not work.

### §1–§2 — the report layer (9 steps, ~24 min)

Previously absent from this page altogether, which meant the page could not be followed to a
finish: everything below depends on the tables these write.

```bash
# NOT named `R` — that would shadow the R binary for the rest of your shell session.
render_rmd() { Rscript -e "rmarkdown::render('$1', quiet=FALSE)"; }

render_rmd analysis/upstream/Isoform_Level_Quantification.Rmd          # 307 s
render_rmd analysis/upstream/NMD_shortread_dge_fullmodel_2026.5.5.Rmd  #  39 s
Rscript    analysis/isopair/isopair_wrapper/01_prepare_data_mashr.R    #  46 s — expression universe
render_rmd analysis/upstream/Isoform_Landscape.Rmd                     #  41 s
render_rmd analysis/upstream/comparison_analysis.Rmd                   #   8 s
render_rmd analysis/upstream/correlation_analysis.Rmd                  #  15 s
render_rmd analysis/upstream/ct_de.Rmd                                 # 314 s
render_rmd analysis/upstream/transcriptional_output.Rmd                # 313 s
render_rmd analysis/upstream/rbp_sr.Rmd                                # 350 s
```

`01_prepare_data_mashr.R` genuinely belongs at that position, interleaved rather than grouped with
the isopair steps below — that is where it ran.

Cell-type scope is asserted, not filtered: the deposited matrices contain exactly the 26
manuscript samples, so `01_prepare_data_mashr.R` verifies that and stops if the input ever
carries anything else, rather than silently dropping it.

### §3–§4 — the isopair chain (4 steps, ~45 min)

```bash
Rscript analysis/isopair/isopair_wrapper/01b_build_isoform_infrastructure.R  # 949 s — structures, CDS, PTC
Rscript analysis/isopair/isopair_wrapper/00b_build_cds_exons.R              #  24 s — CDS exon table
Rscript analysis/isopair/isopair_wrapper/02_build_profiles_mashr.R          # 1350 s — profiles + pairs
Rscript analysis/isopair/isopair_wrapper/03b_rebuild_cache.R                # 378 s — analysis cache
```

`00b` and `03b` were the two steps W97 recorded as missing from this page — which is why claim
4.4.8 could not be rebuilt by following it. The order matters: `00b` runs *between* `01b` and `02`.

### The feature layer (8 steps, ~31 min)

```bash
Rscript analysis/isopair/isopair_wrapper/05k_utr5_all_isoforms.R   # 102 s
Rscript analysis/isopair/isopair_wrapper/05s_orfik_scan.R          # 961 s — the long one
Rscript analysis/isopair/isopair_wrapper/05t_ref_cds_features.R    # 209 s
Rscript analysis/isopair/isopair_wrapper/05u_paralog_annotation.R  #   5 s — queries Ensembl LIVE on a cold cache
Rscript analysis/isopair/isopair_wrapper/05s_b_orfik_scan_extend.R # 490 s
Rscript analysis/isopair/isopair_wrapper/05r_ref_atg_analysis.R    #  49 s
Rscript analysis/isopair/isopair_wrapper/05k_b_utr5_refaug.R       #  37 s
# export_rds.R NEEDS THESE TWO VARIABLES AND WILL NOT RUN WITHOUT THEM. It resolves
# orfik_scan.rds from $NMD_ISOPAIR_CACHE, not from its own directory, so a bare invocation
# stops with "orfik_scan.rds not found" after every step above has succeeded.
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R                               #   9 s — the 8 feature tables
```

`05s` was mis-priced at 75 minutes in `docs/rerun_log.tsv`; that entry records a run that *died*.
It is 961 s.

### The figures the §4 report embeds (10 panels + 2 supplemental, ~4 min)

**THIS SECTION DID NOT EXIST UNTIL 2026-08-11 AND THE §4 REPORT COULD NOT RENDER WITHOUT IT.**
`05_final_report_gencode_scope_2026-07-11.Rmd` calls `include_graphics()` on **14 PNGs**. Every
producer ships, and this page invoked none of them — so the report stopped at chunk 35 of 98 with
*"Cannot find the file(s)"*, after every analysis step above had succeeded. Found by running this
page end to end in a clean room, which is the only way it could have been found: the figures are
already on disk on any machine that has rendered them before.

Each panel script reads the exports written by `data_export.R` in its own directory, so run those
first.

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

# the composites (Figure 3 and Figure 4 themselves)
python3 figures/multipanel/figure3_isopair_and_ptc/figure3_composite.py
python3 figures/multipanel/figure4_ptcneg_and_model/figure4_composite.py
```

**Order matters and the failure is loud, not silent** — a missing PNG stops `knitr` rather than
producing a report with gaps, which is the right behaviour and is why this was caught at all.

### §3 and §4 — the two reports

```bash
render_rmd analysis/upstream/productive_response.Rmd                                       # 855 s — §3
render_rmd analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd #  42 s — §4
```

Both were absent from this page. The §4 report failed twice before it first ran end to end — on a
missing `figures/` output directory, then on figure substrate the projection did not carry (W133) —
so if it fails for you, check those before suspecting the data.

### Supplemental figures for §1–§4

Fourteen producers, run in order by one generator. They read **only** what the chain above
produces plus the deposit — as of 2026-07-30 none of them reads the authoring machine's
`nmd_fig_data/` bundle any more (W139).

```bash
cd analysis/upstream/Figures
Rscript -e 'rmarkdown::render("make_supplemental_figures.Rmd")'
```

**Run it from that directory.** The scripts resolve their inputs through `config/paths.yml` by
walking up to the repository root, but their *output* paths are relative — `supplement_figures/`
beside the generator — so the working directory decides where figures land.

**Its inputs come from two places, and one of them is not yet documented on this page.** The
generator needs `structures.rds` and `cds.rds` from `01b` above, and it needs the DGELists
(`dge_isoform_longread_2026.3.3.rds`, and the `_filtered_` variant) plus the per-cell-type mashr
outputs under `tmp/out/mashr_gene/` and `tmp/out/mashr_isoform/`. Those come from the **§1–§2
upstream reports**, which this page does not list — the same omission W97 records for `00b` and
`03b`. Until it does, run the upstream reports first.

**Two things that will look wrong and are not.** The supplemental producers are `necessary=False`
in `docs/views_G2.tsv` and carry no `run_index`. That is not a gap in the map: `necessary` is a
backward closure from *attributed claims*, and the claim inventory holds no Supplement entries, so
the walk never reaches them. Their code has always been in the graph. And `_nmd_machinery_report.R`
sits in this directory but is **deliberately not run** — it is an internal diagnostic with no
figure number that has never built; see the note at its head.

**Their numbers are not verified against the published supplemental figures.** This step exists so
the supplement rebuilds from the same data as the main text, not to re-check it. The two differ:
the authoring bundle's unfiltered DGEList carried 30,280 isoforms the rebuild does not, and its
filtered one differs in counts across an identical isoform set (W139).

### §5 — the deep-learning model, rebuilt deposit-native

§5 is the one section that does not run from this repository alone. It needs the model
repository — **`NMD_orf_model_v5_4ct`**
([10.5281/zenodo.21536501](https://doi.org/10.5281/zenodo.21536501)) — and a GPU. Clone it
beside this one; `<model repo>` below is its path, and the `python` steps run from inside it.

**You do not need any of this to reproduce §5's reported numbers.** The deposit's `model.zip`
carries the trained checkpoint, the per-isoform predictions and the full interpretability export,
and the §5 figure scripts read those. This chain is for rebuilding them from the deposit.

**Run the wrappers, not bare command lines.** The model repository carries a `slurm_*_dn.sh`
driver for each stage — `_dn` meaning deposit-native — and those are canonical for three reasons a
typed command cannot match: they default `RESULTS_DIR` to the deposit-native tree, they derive the
window from `paths_config.py --selected-tag` rather than repeating a literal, and they carry the
flags whose omission fails silently or late. The command each one runs is shown beneath it so the
page can be read without a cluster, but the wrapper is the thing to submit.

**Three result directories, and two of them are traps.**

| directory | what it is |
|---|---|
| `results_4ct` | the **published** run. Writing here destroys the artifacts the rebuild exists to be compared against, and it is the default of every script that takes `--results-dir` |
| `results_4ct_dn` | **deprecated.** `config_dn.yaml` records that its HDF5 was built from *Channing* inputs — the driver omitted `--config`, so `load_config` fell back to `config.yaml` and `/projects/talisman` resolved on Explorer, so it succeeded — and the tree now sits under `deprecated_2026-08-04/` |
| `results_deposit_h5_2026-08-04` | **the deposit-native tree.** What `config_dn.yaml`'s `hdf5_path` points at and what every `_dn` wrapper defaults to |

**THE WINDOW IS 1000/1000 AND THE CONFIG DEFAULT IS NOT. Read this before scoring anything.** The
deposited checkpoint is `best_model_atg1000_stop1000_seed42.pt`, and the paper agrees with it: §5
reads "1000nt sequence windows around the translational start and stop sites", and reports AUC 0.93
and AUPRC 0.82 on n = 10,522 isoforms with 2,405 NMD susceptible — which is what the deposited
checkpoint scores. `window_size` is a FULL width, not a half-width: `data_prep.py` takes
`half_win = win_size // 2` and reads `half_win` either side of the site, so 1000 means ±500.

`config_dn.yaml` carries the window twice: `window_size_atg: 1000` under `selected:`, and
`window_size_atg: 100` under `data:`, which is the sweep grid's STARTING POINT. `03_train.py` reads
the second — `ws_atg = atg_window or config["data"]["window_size_atg"]` — so a run that does not
pass the window trains at 100/1000, which is neither the selection nor anything the manuscript
reports. The wrappers resolve this by reading the selected tag; anything you write yourself must
pass `--atg-window 1000 --stop-window 1000` explicitly.

**And the model does not take bare sequence** — `data_prep.py`'s nine channels include
exon-junction positions, and the structural features include `n_downstream_ejc`, so scoring an
arbitrary RNA sequence is not a well-posed use of it without transcript structure.

#### The chain

**1 — the feature tables (local, no GPU).** The ORF scan and structural features come from this
repository; `export_rds.R` turns them into the eight TSVs `data_prep.py` reads. These are the same
steps as *The feature layer* above — run them once, and give `export_rds.R` the deposit-native
target rather than its default:

```bash
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R --results-dir results_deposit_h5_2026-08-04
```

**2 — build the HDF5.** `sbatch slurm_build_h5_dn.sh`

```bash
python data_prep.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 --workers 8
```

**3 — patch the stop-codon column.** `sbatch slurm_patch_selected_orfs.sh` — **but read this
first.** `patch_stop_codon.py` fixes an off-by-one in `selected_orfs.tsv`'s `stop_codon` column and
is required for §4.1's χ² test. It is the one step with no deposit-native wrapper: its
`--results-dir` defaults to `results_4ct`, the published tree, and the wrapper passes no arguments
at all. Without `--in-place` it writes `selected_orfs_fixed.tsv` beside the original rather than
replacing it.

```bash
python scripts/patch_stop_codon.py --results-dir results_deposit_h5_2026-08-04 --in-place
```

**4 — the determinism gate.** It must pass before training, and at the real window config.

```bash
python verify_determinism.py --atg 1000 --stop 1000
```

**5 — train.** `sbatch slurm_train_dn.sh` — which also scores `val_clean` afterwards, the
development split.

```bash
TAG=$(python paths_config.py --selected-tag --config config_dn.yaml)   # atg1000_stop1000
python 03_train.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000
python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000 --member-seed 42 --split val_clean
```

**6 — the interpretation layer**, which runs over the full cohort (D74), not the test split.

```bash
sbatch slurm_uorf_dn.sh              # infer_uorf_attention.py, then the uORF metrics
sbatch slurm_deepshap_joint_dn.sh    # deepshap.py --branches joint, 5 runs
sbatch slurm_deepshap_structural_dn.sh
sbatch slurm_kernel_shap_dn.sh       # 11_kernel_shap_branches.py
sbatch slurm_export_chain_dn.sh      # 06_, 07_, 08_, 09_ GC/junction/polyA, 09b_, 09c_, 09d_
```

`slurm_deepshap_all_dn.sh` does joint, structural and `atg stop` in one job if you prefer. The five
DeepSHAP runs vary the background seed — the replicate mean and sd are the uncertainty statement,
which is why there are five and not one.

**7 — the final evaluation. Run this ONCE, and last.**

```bash
sbatch slurm_eval_final_dn.sh
# python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
#     --atg-window 1000 --stop-window 1000 --member-seed 42 --split test_clean --final
```

The test split is untouched through everything above: training early-stops on `val_clean`, the
interpretation chain runs `--split all --full-cohort`, and every export reads the pooled
predictions. `--final` is the affirmation that this is the single evaluation allowed to touch
chr1/3/5/7, and the metrics JSON records `evaluation_class=final_test` so no later reader can
mistake a development number for it.

**8 — the report and the §5 figure data.**

```bash
NMD_RESULTS_DIR=results_deposit_h5_2026-08-04 \
  Rscript -e 'rmarkdown::render("<model repo>/orf_model_report_v5.Rmd",
                                knit_root_dir = normalizePath("<model repo>"))'

# BOTH VARIABLES ARE REQUIRED. MODEL_RESULTS points at the model's own results; FEATURES points
# at what export_rds.R wrote with --results-dir, and the two need not be the same directory.
# SF40's data_export.R needs the same one.
NMD_MODEL_RESULTS=<model repo>/results_deposit_h5_2026-08-04 \
NMD_FEATURES=<model repo>/results_deposit_h5_2026-08-04 \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R
```

Everything below about this chain is load-bearing, and each was learned by getting it wrong:

- **`--member-seed` is required by every consumer of the checkpoint.** `03_train.py` writes
  `best_model_{tag}_seed{N}.pt` — one file per ensemble member — so `evaluate.py`,
  `11_kernel_shap_branches.py` and `deepshap.py` refuse to load an unqualified
  `best_model_{tag}.pt` when members exist, rather than guess which one you meant. Omit it and
  the step dies with `FileNotFoundError` **after** training has succeeded. `--member-seed`
  names the *checkpoint*; `--seed` is the *RNG* seed and is a different thing in the same
  command line.
- **`evaluate.py` requires `--split`, and the test split additionally requires `--final`.** Use
  `--split val_clean` while developing. The test set is scored **once**, deliberately, at the end —
  `--final` is what marks that one evaluation, and the metrics JSON records
  `evaluation_class` so a development number can never be read as a held-out one.
- **Determinism is enforced, not assumed, and it is bounded.** `set_seed` sets
  `cudnn.deterministic`, disables `benchmark`, pins `CUBLAS_WORKSPACE_CONFIG` and turns on
  `torch.use_deterministic_algorithms`. Two seeded runs are then bit-identical *on the same
  GPU model and library versions* — verified on a Tesla V100-SXM2-32GB with torch 2.5.1+cu121.
  A different GPU may still differ, and `NMD_ALLOW_NONDETERMINISM=1` disables all of it
  without leaving a trace in the outputs, so record the environment alongside any run.
- **`--results-dir` is not optional, and it is not only the training scripts.** `03_train.py`,
  `evaluate.py`, `11_kernel_shap_branches.py` and `deepshap.py` all defaulted to a hardcoded
  `results_4ct` and would overwrite the published artifacts the new run exists to be compared
  against. So did `export_rds.R`, `infer_uorf_attention.py`,
  `compute_uorf_attention_metrics.R` and `orf_model_report_v5.Rmd`, which were parameterised
  later — every one of them defaults to the published run, so an invocation that omits the
  override reads or writes the wrong vintage and says nothing about it.
- **`export_rds.R` is the one where omitting it does the most damage,** because it WRITES the
  eight feature tables `data_prep.py` reads. Left at its default it writes into `results_4ct`
  — overwriting the published tables — while `data_prep.py --results-dir results_deposit_h5_2026-08-04`
  reads a directory those tables never reached.
- **The report takes an env var, not a flag.** `rmarkdown::render()` does not pass arguments
  through to the document, so `commandArgs()` inside a knit sees the render call rather than
  the report's own options. `NMD_RESULTS_DIR` is the override.
- **`relabel_tx_summary_4ct.R` is not needed.** It is a proven no-op on a deposit-native
  scaffold (42,043 rows in and out, 0 dropped, 0 flipped), so the dependency on the
  non-deposit `/projects/talisman` mashr root can be dropped without changing a label.
- **The determinism gate must pass first,** and it is not free: `AdaptiveMaxPool1d` had no
  deterministic CUDA backward and was replaced by the mathematically identical
  `x.max(dim=-1)`. Test at the REAL window config — but **not for the reason this line used to
  give.** It said "at ATG=100 `mid_pool` is `Identity`, so the gate passes without exercising the
  op the real run uses." That is **wrong, not merely incomplete** (C69): `ORFEncoder.__init__`
  builds **two** `SequenceCNN` instances, `model.py:102` with `window_size_atg` and `:104` with
  `window_size_stop`, and each sets its own `mid_pool` at `:56` from *its own* window. At
  atg=100/stop=1000 the STOP branch's pool is `MaxPool1d(4)`, so the op *is* exercised. What a
  default-config run leaves untested is the **START branch's pool and the larger kernels** — which
  is still a reason to run at the real config. `verify_determinism.py:109-110` takes
  `--atg`/`--stop` for exactly this.
- **Two things a determinism PASS does not cover,** so do not read one as clearance:
  `verify_pool_equivalence()` is **CPU-only** (`model.py:284` is `torch.randn(n, channels, L)` with
  no device) even though the substitution exists *because*
  `adaptive_max_pool2d_backward_cuda` has no deterministic kernel — it checks the one backend where
  the problem cannot arise. And `verify_determinism.py` contains **no** `autocast`/`GradScaler`
  reference while `03_train.py` trains under AMP, so the real GPU training path is untested. The
  harness itself warns at `:124-125` that a CPU pass proves nothing.
- **DeepSHAP cannot run deterministically at all.** `shap`'s DeepLIFT routes MaxPool
  gradients through `max_unpool1d`, which has no deterministic kernel, so it needs
  `NMD_ALLOW_NONDETERMINISM=1`. That is correct rather than a workaround: the five replicates
  vary the background seed, so the replicate mean and sd — not bitwise determinism — are the
  uncertainty statement.
- **Joint and structural DeepSHAP write the same filename** and silently overwrite or, worse,
  get globbed and averaged together by the consumers. Match the mode of whatever you are
  comparing against; `--branches joint` is what the published Figure 5D used. See W52.

What the deposited run produces, read from the deposit's own
`model/metrics_atg1000_stop1000_seed42_test_clean.json` rather than restated: **AUC 0.9257,
AUPRC 0.8175, n_eval 10,522 with 2,405 NMD susceptible**, `split = test_clean`,
`evaluation_class = final_test`, `best_epoch = 5`, window 1000/1000, `member_seed = 42`. Those
are the numbers §5 of the paper reports at its own precision (0.93 and 0.82). Every §5 claim
whose values move is marked `differs` in `docs/VERIFIED.md`.

`test_clean` is the paralog-free held-out set — chr1/3/5/7 with the 122 paralog-straddling
transcripts carried under their own `test_paralog` label, not filtered out at scoring time. It is
the split the manuscript describes and the published run used, and it is named through
`_split_mask` rather than reconstructed, so `--split test` is a different population and will not
give these numbers.

## What is verified, and what is not

Generated, never asserted:

- [`docs/VERIFIED.md`](docs/VERIFIED.md) — what the rebuild has been observed to do, per claim,
  including which claims carry a value that differs from the published one

## Known gaps

- **THE CLAIM-VALUE TABLE IS WRITTEN TO A DIRECTORY NOTHING COLLECTS.** Ten of the `.Rmd` steps
  below call `claim_emit()`, which writes one row per reported quantity — value, n, population and
  producer — to `$NMD_CLAIM_VALUES`, or to `tmp/claim_values.tsv` when that variable is unset. It
  is unset on the path this page describes, and `tmp/` is git-ignored, so **the file appears, is
  never collected, and is lost with the working directory.** That table is the machine-readable
  record of what your run computed; grepping numbers out of the rendered HTML instead is the method
  that voided 25 outcomes in our own verification. If you want it, set it before you start:

  ```bash
  export NMD_CLAIM_VALUES="$PWD/logs/claim_values.tsv"
  ```

  *Found 2026-08-11 by the paper-validation seat, reading the shipped instructions rather than our
  internal harness — where the same defect had already been fixed, in scripts this package does not
  contain.*

- **YOUR FIGURES WILL NOT LOOK BYTE-IDENTICAL TO THE PAPER'S, AND THAT IS EXPECTED.** The
  published figures were rendered on macOS in **Arial**. The container has no Arial — the font
  is licensed and cannot be redistributed in a public image — so it renders **Liberation Sans**,
  which is a metric-compatible Arial clone. *Advance widths and therefore layout are identical;
  glyph outlines are not.* Panel positions, wrapping and the layout validator all behave the
  same; a pixel-level diff of any figure will differ everywhere.

  This is why figures are verified by comparing their **data exports** rather than their images
  (D120). Where a figure directory ships a `data/` export, the rebuilt export is compared
  byte-for-byte against the tracked one; where none exists, the figure carries a written
  disposition instead of a silent pass.

  The font order is set once in `figures/lib/ggplot_style.py`:
  `Arial → Helvetica Neue → Helvetica → Liberation Sans → Nimbus Sans → DejaVu Sans`. The
  metric clones sit ahead of DejaVu deliberately — before they did, container renders fell
  through to DejaVu, which is **wider**, and SF42 passed the layout validator on the authoring
  laptop while failing it in clean-room job 9049007. If you render on a machine that has Arial
  you will match the paper; the container will not, and neither outcome is a reproduction
  failure.

- **`hexbin` is not a declared dependency** but `render_sr_lr_correlation.R` uses
  `stat_binhex()`. Without it the script still exits 0 and writes a PNG with the layer
  missing. Install it — see [`ENVIRONMENT.md`](ENVIRONMENT.md).
- Stages beyond Isopair stage 1 are wired but not verified end-to-end. Which claims that
  covers is in `docs/VERIFIED.md` under *Not measured* — an unmeasured claim is not an
  unaffected one.
- **§5 has `traced = 0` even though 8 of its 12 claims are checked.** That is not a
  contradiction: the model's producers live in the model repo, so the G2 graph has no path
  for them and cannot yet prove a reader could rebuild them from the deposit. The numbers
  were measured by running the chain above directly. Closing it needs the 06-10 export
  scripts wired into the graph.
- **The model chain needs a GPU and cluster access**, so it is the one part of this document
  that cannot be run from the deposit alone on a laptop. The determinism gate is cheap
  (~30 s) and should be run first every time.
