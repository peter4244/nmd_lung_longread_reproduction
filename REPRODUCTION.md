# Reproducing the analyses

This repository reproduces the analyses in Leshem et al. starting from the Zenodo source-data
record. Raw reads are not required.

## 1. Get the inputs

Everything is in one record and nothing is fetched at run time:
[10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336) — the concept DOI, which always
resolves to the latest version. Cite and follow this one rather than a version DOI, which keeps
pointing at the files of the version that minted it.

**The record is open access.** Download from the Zenodo page; no token or account is needed.
12 data objects, 695 MB compressed, 4.65 GB extracted.

The record's own README owns the layout — follow it to arrange the files into
`<deposit>/source_data`. Two points from it decide whether the rest of this page works:

- **Extracting flat leaves every path pointing into a directory that does not exist.** The wiring
  check in §3 catches this and warns `DOES NOT RESOLVE`. It warns rather than stops, so read its
  output rather than trusting that the command ran clean.
- **On a shared cluster, extract inside a batch allocation.** `sqanti.zip` writes 4.2 GB, and a
  login-node process limit can kill `unzip` part-way; the truncated tree then fails the manifest
  in a way that looks like a corrupt download.

Verify from inside `source_data/` — the manifest's paths are relative to it, not to the record
root:

```bash
cd nmd_deposit_2026/source_data && shasum -a 256 -c ../MANIFEST.sha256
```

All 58 entries should report `OK`. `58 listed files could not be read` means the wrong working
directory rather than a bad download.

## 2. Build the environment

**The container is the validated route.** `nmd_1.3.sif`, in the same record, carries R 4.5.2,
Isopair 0.99.4, edgeR 4.8.2, ORFik 1.30.2, pandoc 3.8.3 and both Python environments.
[`ENVIRONMENT.md`](ENVIRONMENT.md) owns the details: how to fetch or build it, what `nmd_1.3.sif`
is relative to the image the analyses actually ran under, and why three flags are mandatory. The
invocation is repeated once here because every command below assumes it:

```bash
apptainer exec --containall --nv \
  --bind "$REPO":/work --bind "$DEPOSIT":/deposit --pwd /work \
  nmd_1.3.sif <command>
```

**Commands on this page are written for a host shell.** To run one inside the container, wrap it
in the invocation above **and** prefix every `NMD_*` variable with `APPTAINERENV_`:
`NMD_ISOPAIR_CACHE=…` becomes `APPTAINERENV_NMD_ISOPAIR_CACHE=…`, set outside the
`apptainer exec`. `--containall` drops the host environment, so an unprefixed variable arrives
unset and the script stops on a missing file that the previous step spent an hour building
successfully. This applies to every `NMD_*` variable below, `NMD_CLAIM_VALUES` included.

**Verify the image after downloading it — neither manifest covers it.** Zenodo publishes an md5
for every file on the record page; compare it against your copy:

```bash
md5sum nmd_1.3.sif        # md5 -q on macOS
```

Build the environments from source only if you cannot run Apptainer. That gives two of them,
because the work ran in two places — the model on a GPU cluster, the figures and reports on a
laptop.

**R**, from the repository root:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore()'
```

[`renv.lock`](renv.lock) pins 189 packages against **R 4.5.2** and **Bioconductor 3.22**, including
`Isopair`, which is not on CRAN — its record points at GitHub `peter4244/Isopair` at commit
`ea33f74`, so `renv::restore()` fetches it without a separate step.

**Python**, one environment per half:

```bash
conda env create -f environment-model.yml     # the deep-learning model
conda env create -f environment-figures.yml   # figures, reports, pandoc
```

Four things these will not do for you:

- **Compilers.** Several Bioconductor packages build from source on Linux. Without a toolchain
  (`build-essential` or equivalent), `renv::restore()` fails partway with a compiler error rather
  than a missing-package error. macOS needs the Xcode command-line tools.
- **pandoc.** Not an R package, so not in `renv.lock`; it is in `environment-figures.yml`. The
  reports are rmarkdown and will not render without it.
- **A GPU.** `environment-model.yml` installs `torch==2.5.1+cu121`. With no NVIDIA GPU, replace
  the two `pip:` lines with plain `torch==2.5.1`. CUDA changes speed and floating-point
  association, not the model.
- **Byte-identical figures.** See *Known gaps*; fonts are the larger reason.

`hexbin` and `mclust` are referenced by shipped code. `renv::restore()` does not install them; the
container and the `Dockerfile` both do.

**If the conda solve cannot reach its channel, check whether your network blocks
`conda.anaconda.org` before debugging the environment file.** Some institutions block it over
Anaconda's licensing terms, and the error reads like a corrupt package cache. The signature is
that the R half of the build succeeds while the Python half fails, because only the conda step
touches that domain.

## 3. Point the repository at the data

Every path resolves through [`config/paths.yml`](config/paths.yml); nothing is hard-coded.

```bash
ln -s /path/to/nmd_deposit_2026 data_deposit
```

That symlink is deliberately not stored in git: a stored one would carry the original author's
folder layout and then fail on your machine while looking as though it had worked.

Resolution order per key is environment variable `NMD_<KEY>` → `config/paths.yml` → relative to
the repository root. A one-off override is `NMD_DEPOSIT=/elsewhere/source_data Rscript analysis/…`.

Check the wiring before running anything:

```bash
Rscript -e 'source("R/load_config.R"); str(nmd_paths())'
```

**`nmd_paths()` checks the seven deposit-resident keys and warns if any is missing:**

```
Warning: nmd_paths(): 7 of 7 deposit paths DOES NOT RESOLVE -- DEPOSIT=…; SQANTI=…
  The deposit must be arranged as <deposit>/source_data; see REPRODUCTION.md §1.
```

`DEPOSIT`, `SQANTI`, `ISOCALL`, `ANNOT`, `MODEL_RESULTS`, `FEATURES` and `PHENO` must point into
your deposit. **It warns rather than stops**, because rebuilding figures from cached intermediates
is a legitimate run that does not need the deposit mounted — so a warning here is not something to
scroll past. `CACHE`, `OUT` and `ISOPAIR_OUT` are outputs, created on demand and not checked. The `LEGACY`/`FIGDATA`/`MASHR_*`/`ISOPAIR` roots point at
intermediates on the authoring machine; they are not needed to reproduce from the deposit and will
not resolve elsewhere.

## 4. Which count matrix?

The deposit ships two long-read isoform matrices because the published analyses used two isoform
universes. They differ by a single isoform, and **pairing the wrong matrix with the wrong GTF
changes results silently.**

| Section | Counts | Annotation |
|---|---|---|
| §1–2 gene- and isoform-level DE | `nmd_lungcells_counts_4ct.csv` | `sqanti/nmd_lungcells_filtered.gtf` |
| §3–4 isoform-pair / PTC | `nmd_isocall_counts_4ct.csv` | `nmd_isocall_4ct.gtf.gz` |

The shipped code already selects correctly.

## 5. Run

Run the steps in the order below. Outputs go to `tmp/out/` and figures to `figures/`; nothing
writes into the input tree.

Wall times are from one laptop run, so that a twenty-minute step is not mistaken for a hang.
**Expect 2–3× longer on a shared cluster.**

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

**Order matters — `00b` runs *between* `01b` and `02`.**

```bash
Rscript analysis/isopair/isopair_wrapper/01b_build_isoform_infrastructure.R  #  949 s — structures, CDS, PTC
Rscript analysis/isopair/isopair_wrapper/00b_build_cds_exons.R              #   24 s — CDS exon table
Rscript analysis/isopair/isopair_wrapper/02_build_profiles_mashr.R          # 1350 s — profiles + pairs
Rscript analysis/isopair/isopair_wrapper/03b_rebuild_cache.R                #  378 s — analysis cache
```

### 5.3 The feature layer (8 steps, ~31 min)

The last step, `export_rds.R`, lives in the model repository (§5.7) and writes the eight feature
tables the model reads.

```bash
Rscript analysis/isopair/isopair_wrapper/05k_utr5_all_isoforms.R   # 102 s
Rscript analysis/isopair/isopair_wrapper/05s_orfik_scan.R          # 961 s — the long one
Rscript analysis/isopair/isopair_wrapper/05t_ref_cds_features.R    # 209 s
Rscript analysis/isopair/isopair_wrapper/05u_paralog_annotation.R  #   5 s — reads the deposit; no network
Rscript analysis/isopair/isopair_wrapper/05s_b_orfik_scan_extend.R # 490 s
Rscript analysis/isopair/isopair_wrapper/05r_ref_atg_analysis.R    #  49 s
Rscript analysis/isopair/isopair_wrapper/05k_b_utr5_refaug.R       #  37 s

# All three settings are required. export_rds.R resolves orfik_scan.rds from $NMD_ISOPAIR_CACHE,
# not from its own directory, and its --results-dir defaults to results_4ct, the PUBLISHED tree.
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R --results-dir results_deposit_h5_2026-08-04
```

**`export_rds.R` needs a batch allocation, not a login node.** It loads `orfik_scan.rds`
(2.3M rows) into memory and is killed with exit 137 on a login node's per-user cap; 25–30 s under
`--mem=96G`. Every other step in this section runs on a laptop.

### 5.4 The figures the §4 report embeds (~4 min)

The §4 report calls `include_graphics()` on **14 PNGs — not ten — and halts** if one is missing,
so build them first. Each panel script reads the exports written by `data_export.R` in its own
directory.

**5.4 and 5.5 interleave; expect to alternate between them.** Three of the exports these figure
scripts need are written by chunks *inside* the §4 report, and two of those sit after the first
image embed, so the report cannot reach them until the figures exist. Run 5.4, run 5.5 and let the
report halt, run 5.4 again, and repeat: each pass gets further, and the report renders on the pass
where all fourteen PNGs are present. Four passes is typical. **A halt part-way through 5.5 is the
expected intermediate state, not a failure.**

`panel_e_compute.R` is deliberately not listed: it is a fragment that
`figure3…/data_export.R` sources, and running it standalone fails.

```bash
# per-figure data exports
Rscript figures/multipanel/figure3_isopair_and_ptc/data_export.R
Rscript figures/multipanel/figure4_ptcneg_and_model/data_export.R
Rscript figures/supplemental/SF33_CdsAnd3UTR_GENCODE/data_export.R
Rscript figures/supplemental/SF34_TD2Bias_broad/data_export.R

# the ten main-figure panels
for f in figures/multipanel/figure3_isopair_and_ptc/figure3_panel*.py \
         figures/multipanel/figure4_ptcneg_and_model/figure4_panel*.py; do python3 "$f"; done

# the four supplemental panels the report also embeds. SF35 and SF36 ship no data_export.R by
# design: SF35 reads SF34's exports, SF36 reads SF33's.
python3 figures/supplemental/SF33_CdsAnd3UTR_GENCODE/figure_sf33.py
python3 figures/supplemental/SF34_TD2Bias_broad/figure_sf34.py
python3 figures/supplemental/SF35_TD2Bias_occult/figure_sf35.py
python3 figures/supplemental/SF36_CdsAnd3UTR_refAUG/figure_sf36.py

# the composites — Figure 3 and Figure 4 themselves
python3 figures/multipanel/figure3_isopair_and_ptc/figure3_composite.py
python3 figures/multipanel/figure4_ptcneg_and_model/figure4_composite.py
```

### 5.5 The §3 and §4 reports

`render_rmd` is the shell function defined at the top of 5.1; redefine it if you are in a fresh
shell.

```bash
render_rmd analysis/upstream/productive_response.Rmd                                      # 855 s — §3
render_rmd analysis/isopair/isopair_wrapper/05_final_report_gencode_scope_2026-07-11.Rmd  #  42 s — §4
```

### 5.6 Supplemental figures for §1–§4

**SF1–SF23 are not included in this package.** They are produced by
`analysis/upstream/Figures/make_supplemental_figures.Rmd` and its fourteen producers, which are
not part of the reader-facing export.

The supplemental figures that do ship are a different set — the nineteen directories under
`figures/supplemental/`, SF25–SF43. Step 5.4 above builds four of them (SF33–SF36) because the §4
report embeds those. **The remaining fifteen are built in §5.8.**

### 5.7 §5, the deep-learning model

§5 is the one section that does not run from this repository alone. It needs the model repository
**`NMD_orf_model_v5_4ct`** ([10.5281/zenodo.21536501](https://doi.org/10.5281/zenodo.21536501))
and a GPU. Clone it beside this one; `<model repo>` below is its path.

**You do not need to run any of this to reproduce §5's reported numbers.** The deposit's
`model.zip` carries the trained checkpoint, the per-isoform predictions and the full
interpretability export, and the §5 figure scripts read those. This chain rebuilds them.

**Submit the wrappers rather than typing the commands.** The model repository carries a
`slurm_*_dn.sh` driver per stage (`_dn` = deposit-native), which defaults `RESULTS_DIR` to the
deposit-native tree and derives the window from `paths_config.py --selected-tag`. The command each
runs is shown so the page can be read without a cluster.

#### This section needs two environments, and the container is only one of them

**The `slurm_*_dn.sh` wrappers run `$PY <script>` directly on the cluster host; they never invoke
the container.** §2 calls the container the validated route, and it is — for the R side. The
Python steps here do not use it.

**Python (steps 2, 4, 5, 6) — a host interpreter with `torch`, `shap` and `h5py`.** Build it from
**this repository's** `environment-model.yml`, not the model repository's `environment.yml`, which
is a looser specification of the same environment: it pins `pytorch` and `pytorch-cuda` by name
rather than by version.

```bash
conda env create -f environment-model.yml          # the file carries `name: nmd_model`
export PY="$(conda run -n nmd_model which python)"
```

Use `conda`, not `micromamba` — Explorer carries the former and neither of the latter.

**Set `PY` explicitly rather than relying on the default**, which is the authoring account's conda
path. On a shared cluster that path may well be executable by you — group traversal is enough — in
which case the wrapper silently runs someone else's environment and reports success. Confirm what
you got:

```bash
echo "$PY"; "$PY" -c "import sys, torch; print(sys.prefix, torch.__version__)"
```

**The guard cannot catch this for you.** The wrappers test `[ -x "$PY" ]` — executability, not
capability. A bare `python3` on `PATH` passes it, the job starts, and it dies at `import torch`,
which reads like a broken installation rather than a wrong interpreter.

**R (steps 1 and 8) — the container**, with the full invocation from §2.

#### Three result directories, two of them traps

| Directory | What it is |
|---|---|
| `results_4ct` | The **published** run, and the default of every script taking `--results-dir`. Writing here destroys what the rebuild is compared against |
| `results_4ct_dn` | **Deprecated** — its HDF5 was built from Channing inputs |
| `results_deposit_h5_2026-08-04` | **The deposit-native tree**, and what every `_dn` wrapper defaults to |

#### The window is 1000/1000 and the config default is not

The deposited checkpoint is `best_model_atg1000_stop1000_seed42.pt`. `config_dn.yaml` carries the
window twice — `window_size_atg: 1000` under `selected:` and `window_size_atg: 100` under `data:`,
the sweep grid's starting point — and `03_train.py` reads the second. **A run that does not pass
the window trains at 100/1000.** Anything you write yourself must pass
`--atg-window 1000 --stop-window 1000`. `window_size` is a full width: `data_prep.py` takes
`half_win = win_size // 2`, so 1000 means ±500.

The model does not take bare sequence — `data_prep.py`'s nine channels include exon-junction
positions, and the structural features include `n_downstream_ejc` — so scoring an arbitrary RNA
sequence is not a well-posed use of it.

#### The chain

**1 — wire the deposit for the model repository.** `config_dn.yaml` names the deposit as
`../data_deposit/...` **relative to the config file**, so the model repository needs its own
symlink — a *second* one, distinct from §3's. Without it `data_prep.py` stops with
`nmd_lungcells_corrected.fasta not found` while the file sits in the deposit.

```bash
ln -s /path/to/nmd_deposit_2026 data_deposit    # from the directory holding both clones
```

**2 — the feature tables** (batch, no GPU). Already run as the last step of 5.3; repeated because
the rest of this chain depends on it.

**3 — build the HDF5.** `sbatch slurm_build_h5_dn.sh`

```bash
python data_prep.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 --workers 8   # 13 min
```

**4 — the determinism gate.** Must pass before training, at the real window config (~1 min).
`--config` is required; omitting it exits 2 before anything is tested. A pass reports that two
seeded runs are bit-identical in losses and weights.

```bash
python verify_determinism.py --config config_dn.yaml --atg 1000 --stop 1000
```

**5 — train.** `sbatch slurm_train_dn.sh`, which also scores the `val_clean` development split.

```bash
python 03_train.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000            # ~10 s/epoch, early-stops around epoch 19
python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000 --member-seed 42 --split val_clean
```

**6 — the interpretation layer**, which runs over the full cohort, not the test split.

```bash
sbatch slurm_uorf_dn.sh              # infer_uorf_attention.py ONLY — see below
sbatch slurm_deepshap_joint_dn.sh    # deepshap.py --branches joint, 5 runs
sbatch slurm_deepshap_structural_dn.sh
sbatch slurm_deepshap_atgstop_dn.sh  # --branches atg stop → deepshap_all_orfs_summary_*
sbatch slurm_interpret_dn.sh         # 04_interpret_attention.py, 05_interpret_structural.py
sbatch slurm_export_motif_logos_dn.sh  # export_joint_motif_logos.py
sbatch slurm_kernel_shap_dn.sh       # 11_kernel_shap_branches.py
sbatch slurm_export_chain_dn.sh      # 06_, 07_, 08_, 09_ GC/junction/polyA, 09b_, 09c_, 09d_
```

All eight are required. Nothing else writes `deepshap_all_orfs_summary_*`, the four `attention_*`
outputs, or the joint motif logos, so omitting one finishes the step with files missing and no
command having failed.

**`slurm_uorf_dn.sh` runs inference only. The metrics are a second command**, and it must wait for
that job to finish because it consumes what the job writes. Without it
`uorf_attention_metrics.tsv` is never written and step 8 stops on the missing file.

```bash
Rscript compute_uorf_attention_metrics.R --results-dir results_deposit_h5_2026-08-04
```

`slurm_deepshap_all_dn.sh` does joint, structural and `atg stop` in one job — use it *instead of*
those three, not in addition. The five DeepSHAP runs vary the background seed; the replicate mean
and sd are the uncertainty statement.

**7 — the final evaluation. Run once, and last.** `sbatch slurm_eval_final_dn.sh`

```bash
python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000 --member-seed 42 --split test_clean --final
```

The test split is untouched by everything above: training early-stops on `val_clean` and the
interpretation chain runs `--split all --full-cohort`. `--final` marks the single evaluation
allowed to touch chr1/3/5/7, and the metrics JSON records `evaluation_class=final_test`.
`test_clean` is the paralog-free held-out set, with the 122 paralog-straddling transcripts under
their own `test_paralog` label; `--split test` is a different population and gives different
numbers.

**8 — the §5 figure data. This step runs from THIS repository, not the model repository.**

```bash
# Both variables are required. MODEL_RESULTS is the model's own results; FEATURES is what
# export_rds.R wrote with --results-dir. --data-dir is needed only if data_mashr is NOT in-tree:
# nmd_data_dir() defaults to ISOPAIR_OUT, which config/paths.yml declares as tmp/out/data_mashr.
NMD_MODEL_RESULTS=<model repo>/results_deposit_h5_2026-08-04 \
NMD_FEATURES=<model repo>/results_deposit_h5_2026-08-04 \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R \
    --data-dir <path to data_mashr>
```

**`<path to data_mashr>` depends on your layout.** If §3 built it in-tree it is
`tmp/out/data_mashr`. If it was staged *beside* the repository, `--bind "$REPO":/work` does not
reach it — a sibling directory is outside the bound tree — so bind it explicitly and pass the
container-side path:

```bash
apptainer exec --containall --nv --bind "$REPO":/work --bind "$DATA_MASHR":/data_mashr --pwd /work \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R --data-dir /data_mashr
```

Bind the **directory**, never a file inside it — Apptainer refuses a bind destination that does
not exist, where Docker creates one.

**Optional: the model report.** `sbatch slurm_render_dn.sh` renders `orf_model_report_v5.Rmd` from
the model repo (~10 min). Nothing else in §5 depends on it; `data_export_refaug.R` above is what
feeds Figure 5. It renders inside `nmd_1.3.sif`, which is not an implementation detail: bare
`Rscript` is not on the Explorer compute-node PATH, and the report needs `pROC` and `ggseqlogo`,
which `nmd_1.2.sif` lacks and `nmd_1.3.sif` carries. The wrapper also sets `DEPOSIT` to
`<deposit>/source_data/model` rather than `source_data`, because `make_architecture_figure.R` globs
that path directly and refuses to draw rather than guess a window size — one name with two
meanings, disagreeing with `config/paths.yml`.

**A deposit-only reader cannot render the report's Section 1.** It describes the window sweep, and
the deposit carries the deposited model's metrics but not the sweep's 60 files, so that section
fails its guard while every other section renders. Section 9.11 also never renders: it needs a
file no code in either repository produces.

The two sweeps are different runs over different universes, not a current and a stale version, and
their conclusions differ. The published sweep scored all twelve cells on `test_clean`; the
deposit-native one scores `val_clean` across five seeds. The published sweep reported STOP=500 as
optimal by AUC; in the deposit-native sweep the STOP=500 configurations rank 7th, 8th and 9th of
12, and the gap between the best two cells is smaller than the seed-to-seed spread within one cell.

#### Flags whose omission fails silently or late

- **`--member-seed` is required by every consumer of the checkpoint.** `03_train.py` writes
  `best_model_{tag}_seed{N}.pt`, one file per ensemble member, so `evaluate.py`,
  `11_kernel_shap_branches.py` and `deepshap.py` refuse to load an unqualified
  `best_model_{tag}.pt`. Omit it and the step dies with `FileNotFoundError` *after* training has
  succeeded. `--member-seed` names the checkpoint; `--seed` is the RNG seed.
- **`--results-dir` is not optional, and not only for the training scripts.** `export_rds.R`,
  `infer_uorf_attention.py` and `compute_uorf_attention_metrics.R` all default to the published
  run, so omitting it reads or writes the wrong vintage silently.
- **DeepSHAP cannot run deterministically.** `shap`'s DeepLIFT routes MaxPool gradients through
  `max_unpool1d`, which has no deterministic kernel, so it needs `NMD_ALLOW_NONDETERMINISM=1`. The
  five replicates are the uncertainty statement instead.
- **Joint and structural DeepSHAP write the same filename** and will overwrite each other, or be
  globbed and averaged together by the consumers. `--branches joint` is what Figure 5D used.

Determinism is enforced but bounded: `set_seed` pins `cudnn.deterministic`,
`CUBLAS_WORKSPACE_CONFIG` and `torch.use_deterministic_algorithms`, and two seeded runs are then
bit-identical *on the same GPU model and library versions*. A pass does not cover the GPU training
path, which runs under AMP, and `NMD_ALLOW_NONDETERMINISM=1` disables all of it without leaving a
trace in the outputs.

**What the deposited run produces**, read from the deposit's own
`model/metrics_atg1000_stop1000_seed42_test_clean.json`: **AUC 0.9257, AUPRC 0.8175, n_eval 10,522
with 2,405 NMD susceptible**, `evaluation_class = final_test`, `best_epoch = 5`, window 1000/1000,
`member_seed = 42`. These are the numbers §5 reports at its own precision, 0.93 and 0.82.

### 5.8 Figure 5 and the remaining supplemental figures

> **This order is inferred from the scripts, not transcribed from a recorded run.** It was derived
> by reading what each script reads and writes. Treat it as a starting point and check the output,
> particularly for Figure 5, where three different exporters write overlapping panel inputs.

**Every figure directory follows the same shape:** an optional `data_export.R` that writes into a
local `data/` subdirectory, then one figure script that reads it. Script names vary by directory —
`figure_sf37_shap_across_windows.py`, `render_sf25_canonical.R`, and so on — so run what is in
each directory rather than assuming a naming pattern.

**Figure 5** (after §5.7, or from the deposit's `model.zip` alone):

```bash
cd figures/multipanel/figure5_dl_model

python3 data_export_deposit.py          # writes data/panel{B,C,D,E,F,G}_*.tsv

Rscript figure5_panelA_architecture.R   # the architecture diagram; reads no exports
for f in figure5_panel[B-G]*.py; do python3 "$f"; done

python3 figure5_composite.py            # assembles the seven panel PNGs
```

**Panel C reproduces the paper exactly, and the file the panel reads is the right one.** The
exporter writes it twice — `panelC_branch_importance.tsv` over **all isoforms, train and held-out
pooled**, and `panelC_branch_importance_nmd.tsv` over **NMD isoforms only**. The panel script reads
the unsuffixed, pooled file, and that is the population the paper reports: **54.6 / 31.4 / 14.0**
(structural / stop / ATG), which rebuilds from the deposited artifact as 54.64 / 31.41 / 13.95 over
n = 41,776. Use it as the check that your deposit and environment are wired correctly.

The `_nmd` variant gives 59.0 / 29.4 / 11.6 over n = 9,321. It is a legitimate second view, not the
paper's number — do not reconcile the two.

**A caution if you find a different set of percentages elsewhere.** 60.7 / 28.8 / 10.5 belongs to
the **superseded 500 nt model**, NMD-only, n = 2,268; it survives in
`model.published_superseded_2026-07-28/` and in `results_4ct/`. It is not this paper's Figure 5C.

Panel B is held-out test only. Pooling train and held-out is legitimate for interpretation and
never for a performance number.

**Three exporters write Panel G inputs and only one can be right for a given panel.**
`data_export_deposit.py` writes the `_v1` names that `figure5_panelG_uorf_attention.py` actually
reads; `data_export_refaug.R` (§5.7 step 8) writes the `_refaug` variants that SF40 consumes;
`data_export.R` writes a third, unsuffixed set. Confirm which vintage you want before running
more than one — they overwrite into the same `data/` directory.

`make_architecture_figure.R` in this directory is **not** the file of the same name referenced in
§5.7's optional model report; that one lives in the model repository.

**The fifteen remaining supplemental figures.** SF25–SF32 read isopair outputs, so run them after
§5.5. SF37–SF43 read model outputs, so run them after §5.7 or from `model.zip`.

```bash
# after §5.5 — isopair-derived
for d in figures/supplemental/SF2[5-9]_* figures/supplemental/SF3[0-2]_*; do
  [ -f "$d/data_export.R" ] && Rscript "$d/data_export.R"
done   # then run the single figure script in each directory

# after §5.7 — model-derived
for d in figures/supplemental/SF3[7-9]_* figures/supplemental/SF4[0-3]_*; do
  [ -f "$d/data_export.R" ] && Rscript "$d/data_export.R"
done   # then run the single figure script in each directory
```

Of these fifteen, six ship a `data_export.R` — SF27, SF28, SF29, SF30, SF31 and SF40; the rest
read a sibling's exports or the deposit directly.

**SF43 is the exception: it has no `data_export.R` and needs a chain run first.**
`figure_s_model_comparison.py` reads `metrics_summary_<date>.tsv` and
`per_isoform_scores_<date>.tsv`, which nothing in `figures/` produces — they come from
`analysis/predictor_comparison/`, run in order:

```bash
Rscript analysis/predictor_comparison/01_extract_our_isoforms.R
Rscript analysis/predictor_comparison/02_score_nmdetective_b.R
Rscript analysis/predictor_comparison/03_score_nmdep_rule_baseline.R
Rscript analysis/predictor_comparison/04_compute_metrics.R
```

**Mind the datestamp.** `04_compute_metrics.R` names its outputs from a `DATESTAMP` variable, while
`figure_s_model_comparison.py` and `nmd_predictor_comparison.Rmd` each hardcode a different date.
If SF43 reports a missing input, that mismatch is the first thing to check.

## What is verified, and what is not

[`docs/VERIFIED.md`](docs/VERIFIED.md) records, per claim, what the rebuild has been observed to
do, including which claims carry a value that differs from the published one.

## Known gaps

- **Set `NMD_CLAIM_VALUES` before you start, or lose the claim table.** Ten of the `.Rmd` steps
  call `claim_emit()`, which writes one row per reported quantity — value, n, population,
  producer — to `$NMD_CLAIM_VALUES`, or to `tmp/claim_values.tsv` when unset. `tmp/` is
  git-ignored, so the file is written and then lost with the working directory. It is the
  machine-readable record of what your run computed.

  ```bash
  export NMD_CLAIM_VALUES="$PWD/logs/claim_values.tsv"
  # inside the container, export is not enough:
  APPTAINERENV_NMD_CLAIM_VALUES=$PWD/logs/claim_values.tsv apptainer exec --containall …
  ```

- **Your figures will not be byte-identical to the paper's.** The published figures were rendered
  on macOS in **Arial**. The container has no Arial — it is licensed and cannot be redistributed —
  so it renders **Liberation Sans**, a metric-compatible clone. Advance widths and therefore layout
  are identical; glyph outlines are not, so a pixel diff differs everywhere while panel positions
  and wrapping match. Figures are therefore verified by comparing their **data exports**, not their
  images. The font order is set once in `figures/lib/ggplot_style.py`:
  `Arial → Helvetica Neue → Helvetica → Liberation Sans → Nimbus Sans → DejaVu Sans`.

- **Stages beyond Isopair stage 1 are wired but not verified end-to-end.** Which claims that covers
  is in `docs/VERIFIED.md` under *Not measured*; an unmeasured claim is not an unaffected one.

- **The §5 chain needs a GPU and cluster access**, so it is the one part of this document that
  cannot be run from the deposit alone on a laptop.
