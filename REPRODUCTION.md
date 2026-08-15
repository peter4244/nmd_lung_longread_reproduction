# Reproducing the analyses

This repository reproduces the analyses in Leshem et al. starting from the Zenodo source-data
record. Raw reads are not required.

## 1. Get the inputs

| What | Where |
|---|---|
| Source data (12 data objects — 7 files and 5 archives; 695 MB compressed, 4.65 GB extracted) | Zenodo **`10.5281/zenodo.21544336`** — the concept DOI, which always resolves to the latest version. Follow this one, not a version DOI: a version DOI keeps pointing at the files of the version that minted it |
| Tan et al. supplementary tables (§2 only) | already in the record, `tan_2025_supplementary/` |
| MSigDB hallmark gene sets (§3 only) | already in the record, `annotation/msigdb_hallmark_H.tsv` |

Every input is in the record and nothing is fetched at run time.

Zenodo presents a flat file list; the code expects the data one level down. Fetch and assemble:

`10.5281/zenodo.21544336` is the **concept** DOI. Its id is not a record id: Zenodo's
`/api/records/<id>/files` endpoint serves only *version* ids and returns 404 for a concept id, with
or without a token. Resolve the concept to its latest version first:

```bash
# Resolve the concept DOI to its latest VERSION record id. -L is required: the concept id answers
# 302 with an HTML redirect body, which json.load then chokes on.
REC=$(curl -sL https://zenodo.org/api/records/21544336/versions/latest \
      | python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])')
echo "resolved to version record $REC"   # 21878996 at the time of writing

# the record's files are access-restricted; without a token this returns 403
for f in $(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
             "https://zenodo.org/api/records/$REC/files" \
           | python3 -c 'import json,sys; [print(e["key"]) for e in json.load(sys.stdin)["entries"]]'); do
  curl -sL -H "Authorization: Bearer $ZENODO_TOKEN" -o "$f" \
    "https://zenodo.org/api/records/$REC/files/$f/content"
done

mkdir -p nmd_deposit_2026/source_data
mv *.csv *.gtf.gz *.zip nmd_deposit_2026/source_data/
mv MANIFEST.sha256 UPLOAD_MANIFEST.sha256 nmd_deposit_2026/
(cd nmd_deposit_2026/source_data && for z in *.zip; do unzip -q "$z" && rm "$z"; done)
```

Pass the token in the header, never in the URL. If you have no token, the file listing returns
`{"status": 403, "message": "Permission denied."}` and nothing downloads — the record's files are
restricted, which is a permissions problem and not something retrying fixes.

**On a shared cluster, extract in a batch allocation.** `sqanti.zip` writes 4.0 GB, and a login-node
process limit can kill `unzip` part-way; the truncated tree then fails the manifest in a way that
looks like a corrupt download.

Then verify **from inside `source_data/`** — the manifest's paths are relative to it, not to the
record root:

```bash
cd nmd_deposit_2026/source_data && shasum -a 256 -c ../MANIFEST.sha256
```

All 58 entries should report `OK`. If instead you see `58 listed files could not be read`, that is
the wrong working directory rather than a bad download — the manifest's paths are relative to
`source_data/`.

## 2. Build the environment

**The container is the validated route.** `nmd_1.2.sif` in the same Zenodo record carries R 4.5.2,
Isopair 0.99.4, edgeR 4.8.2, ORFik 1.30.2, pandoc 3.8.3 and both Python environments, and is what
every clean-room run of this chain has used. Build the environments from source only if you cannot
run apptainer. See [`ENVIRONMENT.md`](ENVIRONMENT.md).

**Two flags decide whether the container works, and omitting either fails in a way that does not
name it:**

```bash
apptainer exec --containall --nv \
  --bind "$REPO":/work --bind "$DEPOSIT":/deposit --pwd /work \
  nmd_1.2.sif <command>
```

- **`--pwd`** — `--containall` starts you in an empty `$HOME`, not the host working directory, and
  every command on this page is repo-relative. Without it the first step dies with
  `The file 'analysis/upstream/ct_de.Rmd' does not exist` while the file is plainly there.
- **`--nv`** — without it the container sees no GPU, `torch.cuda.is_available()` returns False, and
  §5.7 **trains on CPU without an error**, emitting a plausible checkpoint by a different numerical
  path under mixed precision. This one is silent, which makes it the more dangerous of the two.

Zenodo publishes an md5 for the image; neither manifest covers it, so check it explicitly:

```bash
curl -s https://zenodo.org/api/records/21544336/files \
  | python3 -c 'import json,sys; print([e["checksum"] for e in json.load(sys.stdin)["entries"] if e["key"]=="nmd_1.2.sif"][0])'
md5sum nmd_1.2.sif
```

Building from source instead gives two environments, because the work ran in two places: the model
on a GPU cluster, the figures and reports on a laptop.

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

Two packages §5 needs were missing from these files until 2026-08-12 and are now included: `shap`,
without which `deepshap.py` cannot run and Figure 5D cannot be produced, and `r-proc`, without
which the §5 model report cannot render. If you built an environment from an earlier release, add
them.

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

`hexbin` and `mclust` are referenced by shipped code. `renv::restore()` does not install them; the
container and the `Dockerfile` both do.

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
writes into the input tree.

Wall times are from one laptop run and are there so a twenty-minute step is not mistaken for a
hang. **Expect 2–3× longer on a shared cluster** — measured 2.3× in aggregate and 3.3× at worst on
Explorer.

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
(2.3M rows) into memory and is killed with exit 137 on a login node's per-user cap. Measured:
25–30 s under `--mem=96G`. Every other step in this section runs on a laptop.

**Inside the container these two variables need the `APPTAINERENV_` prefix.** `--containall` drops
the host environment, so the plain shell prefix above arrives unset and the script stops with
`orfik_scan.rds not found` — pointing at your data, immediately after the step above spent an hour
building that very file successfully:

```bash
APPTAINERENV_NMD_ISOPAIR_CACHE=/work/tmp/out/data_mashr/analysis_cache \
APPTAINERENV_NMD_ISOPAIR_DATA=/work/tmp/out/data_mashr \
  apptainer exec --containall … Rscript /model/export_rds.R --results-dir results_deposit_h5_2026-08-04
```

The same applies to every `NMD_*` variable this page sets, `NMD_CLAIM_VALUES` and
`NMD_RESULTS_DIR` included.

### 5.4 The figures the §4 report embeds (~4 min)

The §4 report calls `include_graphics()` on 14 PNGs and **halts** if one is missing, so build them
first. Each panel script reads the exports written by `data_export.R` in its own directory.

**5.4 and 5.5 interleave — expect to alternate between them.** Three of the exports these figure
scripts need are written by chunks *inside* the §4 report, and two of those chunks sit after the
first image embed at line 1257, so the report cannot reach them until the figures exist. Neither
section completes in one pass. Run 5.4, run 5.5 and let the report halt, run 5.4 again, and repeat:
each pass gets further, and the report renders on the pass where all fourteen PNGs are present. Four
passes is typical. A halt part-way through 5.5 is the expected intermediate state, not a failure.

`panel_e_compute.R` is deliberately **not** listed below: it is a fragment that
`figure3…/data_export.R` sources, and running it as a standalone step fails.

```bash
# per-figure data exports
Rscript figures/multipanel/figure3_isopair_and_ptc/data_export.R
Rscript figures/multipanel/figure4_ptcneg_and_model/data_export.R
Rscript figures/supplemental/SF33_CdsAnd3UTR_GENCODE/data_export.R
Rscript figures/supplemental/SF34_TD2Bias_broad/data_export.R

# the ten main-figure panels
for f in figures/multipanel/figure3_isopair_and_ptc/figure3_panel*.py \
         figures/multipanel/figure4_ptcneg_and_model/figure4_panel*.py; do python3 "$f"; done

# THE FOUR SUPPLEMENTAL PANELS THE REPORT ALSO EMBEDS — fourteen in total, not ten. Omit these and
# the report halts in the next section, on a figure this page never told you to build. SF35 and SF36
# ship no data_export.R by design: SF35 reads SF34's exports, SF36 reads SF33's.
python3 figures/supplemental/SF33_CdsAnd3UTR_GENCODE/figure_sf33.py
python3 figures/supplemental/SF34_TD2Bias_broad/figure_sf34.py
python3 figures/supplemental/SF35_TD2Bias_occult/figure_sf35.py
python3 figures/supplemental/SF36_CdsAnd3UTR_refAUG/figure_sf36.py

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

**Not included in this package.** These are produced by
`analysis/upstream/Figures/make_supplemental_figures.Rmd` and its fourteen producers, which cover
the SF1–SF23 range and are not part of the reader-facing export. Running the chain above does not
build them, and this page previously gave commands for them that could not work here.

The supplemental figures that **do** ship are a different set — `figures/supplemental/SF25`…`SF43`,
each with its own producer, built as part of step 5.4 and by the §5 chain below.

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

#### The environment for this section is NOT the container, and §2 does not tell you that

**The `slurm_*_dn.sh` wrappers run `$PY <script>` directly on the cluster host. They never invoke
the container.** §2 calls the container "the validated route" and it is — for the R side. It is
not what the Python steps here use, and nothing on this page said so until 2026-08-12. A reader
who followed §2, obtained `nmd_1.2.sif`, and then submitted a wrapper was running neither the
image they had just built nor anything the page had described.

So this section needs **two** environments, and you must supply the first one yourself:

**Python (steps 2, 4, 5, 6) — a host interpreter with `torch`, `shap` and `h5py`.** Build it from
**this repository's** `environment-model.yml`:

```bash
conda env create -f environment-model.yml          # the file carries `name: nmd_model`
export PY="$(conda run -n nmd_model which python)"
```

**`conda`, not `micromamba`.** Explorer has `conda` on `PATH` and has neither `micromamba` nor
`mamba` — checked 2026-08-12. `mamba` is much faster where it exists, but an instruction naming a
tool the target cluster does not carry fails at the reader's first command.

**If the solve cannot reach the channel, check whether your network blocks Anaconda before you
debug the environment file.** Both files use `conda-forge`, which is hosted at
`conda.anaconda.org`, and a number of institutions block that domain over Anaconda's commercial
licensing terms. On a blocked network micromamba retries eight times and then reports
`Download error (35) SSL connect error` and `Subdir conda-forge/noarch not loaded!`, which reads
like a corrupt package cache — its own message suggests `mamba clean -a` — and is nothing of the
kind. Measured on 2026-08-12: `conda.anaconda.org` and `repo.anaconda.com` both timed out on one
institutional VPN while `p3m.dev`, `github.com` and `pypi.org` all answered normally, and all four
answered off it. That is what makes it diagnosable: **the R half of this build will succeed while
the Python half fails**, because the R packages come from P3M and only the conda step touches
Anaconda.

**Use that file and not the model repository's `environment.yml`**, which is a different and looser
specification of the same environment — it pins `pytorch` and `pytorch-cuda=11.8` by name rather
than by version. `environment-model.yml` was read off the Explorer environment that actually
produced the §5 numbers and pins every version. Two files describing one environment is a trap we
have not yet closed; until we do, this is the one of record.

`PY` is the documented override and every `_dn` wrapper honours it. **Set it explicitly rather
than relying on the default**, which is the authoring account's conda path
(`/home/p.castaldi/.conda/envs/nmd_model/bin/python`). On a shared cluster that path may well be
*executable by you* — group traversal is enough — in which case the wrapper silently runs someone
else's environment and reports success. Confirm what you actually got:

```bash
echo "$PY"; "$PY" -c "import sys, torch; print(sys.prefix, torch.__version__)"
```

**And note the guard cannot catch this for you.** The wrappers test `[ -x "$PY" ]` —
*executability, not capability*. A bare `python3` on `PATH` passes it, the job starts, and it dies
at `import torch`, which reads like a broken installation rather than a wrong interpreter.

**R (steps 1 and 8) — the container.** These need R 4.5.2 with `Isopair`, `ORFik`, `rmarkdown` and
pandoc, which is what the image is for. Run them through it, with the full invocation from §2:

```bash
apptainer exec --containall --nv --bind "$REPO":/work --bind "$DEPOSIT":/deposit --pwd /work \
  nmd_1.2.sif Rscript <script>
```

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

**0 — wire the deposit for the model repository.** `config_dn.yaml` names the deposit as
`../data_deposit/...` **relative to the config file**, so the model repository needs its own
symlink beside it — this is a *second* one, distinct from the `data_deposit` created inside this
repository in §3. Without it `data_prep.py` stops with `nmd_lungcells_corrected.fasta not found`
while the file sits in the deposit.

```bash
# from the directory holding both clones
ln -s /path/to/nmd_deposit_2026 data_deposit
```

**1 — the feature tables** (batch, no GPU — see §5.3 on memory). Already run as the last step of
5.3; repeated here because the rest of this chain depends on it.

```bash
NMD_ISOPAIR_CACHE=tmp/out/data_mashr/analysis_cache \
NMD_ISOPAIR_DATA=tmp/out/data_mashr \
  Rscript <model repo>/export_rds.R --results-dir results_deposit_h5_2026-08-04
```

**2 — build the HDF5.** `sbatch slurm_build_h5_dn.sh`

```bash
python data_prep.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 --workers 8   # 13 min
```

**3 — the stop-codon patch. SKIP THIS: it does not apply to the current pipeline.** The model
repository's build order lists `scripts/patch_stop_codon.py` as required for §4.1's χ² test. It
reads `h5["w500/stop_windows"]`. A current build **does** contain a `w500` group, alongside `w100`,
`w1000` and `w2000`, but its datasets are named `atg_codes`/`stop_codes` — the group is fine, the
dataset name is from an earlier schema — so it exits 1 with `KeyError: 'stop_windows' doesn't exist`. It also has nothing to do:
run against a deposit-native build it reports `Pre-fix canonical (TGA/TAA/TAG) stop rate: 100.0%`
before failing. Measured 2026-08-11 against `results_deposit_h5_2026-08-04`.

**4 — the determinism gate.** Must pass before training, at the real window config (~1 min).

```bash
python verify_determinism.py --config config_dn.yaml --atg 1000 --stop 1000
```

`--config` is required; omitting it exits 2 with `the following arguments are required: --config`
before anything is tested.
A pass reports *two seeded runs are bit-identical in losses and weights*.

**5 — train.** `sbatch slurm_train_dn.sh`, which also scores the `val_clean` development split.

```bash
python 03_train.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000            # ~10 s/epoch, early-stops around epoch 19
python evaluate.py --config config_dn.yaml --results-dir results_deposit_h5_2026-08-04 \
    --atg-window 1000 --stop-window 1000 --member-seed 42 --split val_clean
```

**6 — the interpretation layer**, which runs over the full cohort, not the test split.

```bash
sbatch slurm_uorf_dn.sh              # infer_uorf_attention.py ONLY — see the follow-up below
sbatch slurm_deepshap_joint_dn.sh    # deepshap.py --branches joint, 5 runs
sbatch slurm_deepshap_structural_dn.sh
sbatch slurm_deepshap_atgstop_dn.sh  # --branches atg stop → deepshap_all_orfs_summary_*
sbatch slurm_interpret_dn.sh         # 04_interpret_attention.py, 05_interpret_structural.py → attention_*
sbatch slurm_export_motif_logos_dn.sh  # export_joint_motif_logos.py → motif_logo_{atg,stop}_joint_*
sbatch slurm_kernel_shap_dn.sh       # 11_kernel_shap_branches.py
sbatch slurm_export_chain_dn.sh      # 06_, 07_, 08_, 09_ GC/junction/polyA, 09b_, 09c_, 09d_
```

**Three of those were missing from this list until 2026-08-14, and they are exactly the producers
of the seven interpretation files that then look "absent".** `slurm_deepshap_atgstop_dn.sh`,
`slurm_interpret_dn.sh` and `slurm_export_motif_logos_dn.sh` are not optional extras: nothing else
writes `deepshap_all_orfs_summary_*`, the four `attention_*` outputs, or the joint motif logos. A
reader who ran only the six previously listed here would finish the step with those files missing
and no command having failed.

**`slurm_uorf_dn.sh` RUNS INFERENCE ONLY. The metrics are a second command**, and it must wait for
that job to finish — it consumes what the job writes:

```bash
# AFTER slurm_uorf_dn.sh completes. --results-dir is required: the script defaults to the
# published run, so omitting it silently reads the wrong vintage (see the flags note below).
Rscript compute_uorf_attention_metrics.R --results-dir results_deposit_h5_2026-08-04
```

Without it `uorf_attention_metrics.tsv` is never written, and step 8 stops on the missing file.
`slurm_uorf_dn.sh` names this script in a comment at its foot but has no runnable invocation of it,
which is why the gap was invisible: every command in the step reported success.

`slurm_deepshap_all_dn.sh` does joint, structural and `atg stop` in one job — use it *instead of*
those three, not in addition. The five DeepSHAP runs vary the background seed; the replicate mean
and sd are the uncertainty statement.

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

**8 — the §5 figure data.**

**THIS STEP RUNS FROM THIS REPOSITORY, NOT THE MODEL REPOSITORY.** Steps 2–7 are `sbatch` from
the model repo; this one is an `Rscript` from the reproduction tree. The paths below are the only
thing that still points at the model repo.

```bash
# Both variables are required. MODEL_RESULTS is the model's own results; FEATURES is what
# export_rds.R wrote with --results-dir. SF40's data_export.R needs the same one.
# --data-dir IS REQUIRED IN A CLEAN ROOM. nmd_data_dir() falls back to .P$OUT/data_mashr, which
# exists in the authoring tree and does not exist in a fresh checkout; it fails loudly rather
# than falling back to another vintage, so the step stops here without it.
NMD_MODEL_RESULTS=<model repo>/results_deposit_h5_2026-08-04 \
NMD_FEATURES=<model repo>/results_deposit_h5_2026-08-04 \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R \
    --data-dir <path to data_mashr>
```

**`<path to data_mashr>` DEPENDS ON YOUR LAYOUT, AND THE TWO LAYOUTS DIFFER IN A WAY THAT BITES.**
If §3 built it in-tree it is `tmp/out/data_mashr`. If it was staged *beside* the repository rather
than inside it — which is how the clean room handed to an outside account was arranged — then it is
an absolute path, **and `--bind "$REPO":/work` does not reach it**: a sibling directory is outside
the bound tree, so inside the container the path resolves to nothing. Bind it explicitly and pass
the container-side path:

```bash
apptainer exec --containall --bind "$REPO":/work --bind "$DATA_MASHR":/data_mashr --pwd /work … \
  Rscript figures/multipanel/figure5_dl_model/data_export_refaug.R --data-dir /data_mashr
```

Bind the **directory**, never a file inside it — Apptainer refuses a bind destination that does not
exist, where Docker creates one, so a file-level bind that works locally fails in the `.sif`.

**`orf_model_report_v5.Rmd` NOW RENDERS, AND ON 2026-08-15 IT DID SO END TO END FOR THE FIRST
TIME.** It is an optional step — nothing else in §5 depends on it — but it is no longer a step
you *cannot* run, which is what this page said for two days. The history is kept because the old
reason is still quoted in older renders.

```bash
sbatch slurm_render_dn.sh        # from the MODEL repo, ~10 min, writes render_out/
```

**IT RENDERS INSIDE `nmd_1.3.sif`, AND THAT IS NOT AN IMPLEMENTATION DETAIL.** Bare `Rscript` is
not on the Explorer compute-node PATH — a plain `Rscript -e rmarkdown::render(...)` dies with
`command not found`, exit 127 — and the node's module R would not help either, because this
report needs `pROC` and `ggseqlogo`, which `nmd_1.2.sif` lacks and `nmd_1.3.sif` carries. The
wrapper also redirects `TMPDIR` off the 64 MB `--containall` tmpfs, and sets `DEPOSIT` to
`<deposit>/source_data/**model**` rather than `source_data`, because `make_architecture_figure.R`
globs that path directly and **refuses to draw rather than guess a window size** if it resolves
nothing. That last one disagrees with `config/paths.yml`, which defines `DEPOSIT` one level up —
one name, two meanings.

**And the disagreement is masked for us and not for a reader.** The figure script tries three
roots and the first that resolves wins: `$DEPOSIT`, then `./data_deposit/source_data/model`, then
a hardcoded home path. Run from the repo root the second root resolves, so the wrong `$DEPOSIT`
never surfaces here — which is presumably why nobody noticed. It fails only when roots 2 and 3
also miss: a reader who extracted the deposit elsewhere and set `DEPOSIT` per `paths.yml`, which
is exactly the scenario the variable exists for.

**SECTION 9.11 NEVER RENDERS, BY CONSTRUCTION.** It needs
`deepshap_all_orfs_summary_<tag>_seed42.tsv`, which **no code in either repository produces** —
the only mention is a comment — and which exists in no results tree. The report guards it with
`file.exists`, so the section reports itself as intentionally not generated. Every other input the
report needs now exists.

*What was true until 2026-08-14, and is now fixed:* the report read `metrics_atg{A}_stop{S}.json`,
one per cell of the **published** 2026-04 sweep. This chain names its outputs
`metrics_<tag>_seed<N>_<split>.json`, so none of the twelve existed, `metrics_list` came back
empty, and the render died on `object 'auprc' not found` — a message naming neither the absent
files nor the reason. Measured 2026-08-13, job 9162085.

*What is true now:* on Pete's instruction the report's Section 1 reads the **deposit-native**
sweep, `results_sweep_dn_2026-08-04` — 12 configurations × 5 paired seeds, scored on `val_clean`.
It resolves that tree through `NMD_SWEEP_DIR`, separately from `NMD_RESULTS_DIR`, because the
sweep and the selected model are different runs in different directories.

**THE ONE THING A DEPOSIT-ONLY READER STILL CANNOT DO IS SECTION 1.** The deposit carries two
metrics files — the superseded `metrics_atg500_stop500.json` and the deposited model's
`metrics_atg1000_stop1000_seed42_test_clean.json` — and **not** the sweep's 60. Section 1
describes the window sweep, so on a deposit-only tree it fails its guard while every other
section renders. Running the sweep to recover it is a far larger job than the whole of §5.7; the
60 files are small and depositing them would close it.

**The two sweeps are different runs over different universes, not a current version and a stale
one — and their conclusions differ.** The published one scored all twelve cells on `test_clean`,
which is the leak D67 records; the deposit-native one scores `val_clean` across five seeds. The
published sweep reported *"STOP=500 is optimal by AUC"*; in the deposit-native sweep the STOP=500
configurations rank 7th, 8th and 9th of 12 and the gap between the best two cells is smaller than
the seed-to-seed spread within one cell. **Do not reconcile the two by editing numbers.**

Nothing else in §5 depends on this report: `data_export_refaug.R` above is what feeds Figure 5.

**Flags whose omission fails silently or late:**

- **`--member-seed` is required by every consumer of the checkpoint.** `03_train.py` writes
  `best_model_{tag}_seed{N}.pt`, one file per ensemble member, so `evaluate.py`,
  `11_kernel_shap_branches.py` and `deepshap.py` refuse to load an unqualified
  `best_model_{tag}.pt`. Omit it and the step dies with `FileNotFoundError` *after* training has
  succeeded. `--member-seed` names the checkpoint; `--seed` is the RNG seed.
- **`--results-dir` is not optional, and not only for the training scripts.** `export_rds.R`,
  `infer_uorf_attention.py` and `compute_uorf_attention_metrics.R` all default to the published
  run, so omitting it reads or writes the wrong vintage silently.
- **`orf_model_report_v5.Rmd` defaults to the published run too, and for it that default is
  correct** — it describes the initial window sweep, which only the published tree holds. It
  takes `NMD_RESULTS_DIR` rather than a flag, because `rmarkdown::render()` does not pass
  arguments through to the document. It is not part of this chain; see step 8.
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

  **Inside the container, `export` is not enough.** `--containall` and `--cleanenv` drop the host
  environment, so the variable arrives unset and the table is lost exactly as described above. Pass
  it through explicitly:

  ```bash
  APPTAINERENV_NMD_CLAIM_VALUES=$PWD/logs/claim_values.tsv apptainer exec --containall …
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
