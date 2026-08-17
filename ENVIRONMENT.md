# Environment

The software versions the manuscript analyses ran under, and how to run the container correctly.

## Getting the environment

**Fetch the image (recommended).** The Zenodo record
([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) carries `nmd_1.3.sif`
(8.1 GB) for Apptainer, and the same environment as a Docker archive (8.8 GB) for Docker users:

```bash
docker load -i nmd_docker_1.0_2026-08-13_verified.tar
```

`nmd_1.3.sif` is a **corrected rebuild**, not a byte-identical copy of the image the results came
out of. The analyses ran under `nmd_1.2.sif`, which lacked `pROC` and `ggseqlogo` and whose bare
`Rscript` resolved to a conda R holding none of the analysis packages. 1.3 is rebuilt from the
same pinned lockfile with both fixed, so it is the environment the paper describes.

**`nmd_1.2.sif` is not archived and is not meant to be.** This document previously said it remained
retrievable from an earlier version of the Zenodo record; checked against all six versions on
2026-08-16, no version has ever carried it. Nor should one: it is the image with the missing
packages described above, and offering it beside the corrected build would hand readers an
environment that cannot render the model report. 1.3 is the reproducible environment. What 1.2 was
is recorded in `BUILD_NOTES.md` — an Apptainer layer over 1.1 adding four Liberation Sans fonts,
with no Docker path to it — which is why figures are verified by comparing data exports rather than
by comparing pixels.

**Or build it.** `Dockerfile` in the repository root restores all 189 packages in `renv.lock` and
both conda environments, failing loudly if it cannot. The model environment is included by
default.

```bash
docker build --platform linux/amd64 -t nmd:1.3 .
docker save nmd:1.3 -o nmd.tar
apptainer build nmd_1.3.sif docker-archive://nmd.tar
```

Save a *locally built* image. `docker save` of a **pulled** multi-arch image produces a
multi-manifest tarball that Apptainer refuses; `--provenance` and SBOM attestations reintroduce
the same problem.

## Running the container

The invocation is part of the environment. Two flags are mandatory everywhere, and neither announces
itself when omitted — each produces a plausible wrong answer rather than an error. A third, `--nv`,
belongs only to §5.7.

```bash
export REPO=/path/to/nmd_lung_longread_reproduction   # your clone of this repository
export DEPOSIT_ROOT=/path/to/nmd_deposit_2026         # the deposit root — the directory holding source_data/

export MODEL_REPO=/path/to/NMD_orf_model_v5_4ct   # only if you are running §5.3 or §5.7

apptainer exec --containall \
  --bind "$REPO":/work --bind "$DEPOSIT_ROOT":"$DEPOSIT_ROOT" \
  --bind "$MODEL_REPO":"$MODEL_REPO" --pwd /work \
  nmd_1.3.sif <command>
```

**Bind the model repository too, and bind it at its own host path.** §5.7 tells you to clone it
*beside* this repository — so it is outside `/work` and `--containall` cannot see it. §5.3's last
step runs `Rscript <model repo>/export_rds.R` and is where this bites first, several sections before
§5.7 introduces the repo. Same rule as the deposit: bind it where it already lives, so paths resolve
identically inside and out.

**Add `--nv` only for §5.7.** It is the GPU flag and §5.7 is the only section that uses a GPU. On a
CPU node it is harmless but prints `WARNING: Could not find any nv files on this host!` on *every*
container call — measured across a full §5.1–§5.5 run, that is hundreds of identical warnings in
the logs a reader is being told to read, which is how a real warning gets missed.

**Bind the deposit at its own host path, not at `/deposit`.** REPRODUCTION.md §3 has you create a
`data_deposit` symlink inside the repository pointing at the deposit root. That symlink stores an
absolute host path, so if the deposit is bound anywhere else, `--containall` leaves it dangling
inside the container and all seven deposit keys fail to resolve — the wiring check reports
`DOES NOT RESOLVE` and nothing else works. Binding it at the same path it has on the host makes the
symlink resolve identically inside and out.

**Do not export a variable called `DEPOSIT`.** Two §5 scripts —
`analysis/section5/derive_section5_numbers.py` and
`figures/multipanel/figure5_dl_model/data_export_deposit.py` — read `DEPOSIT` as the **model
subdirectory** (`<deposit>/source_data/model`), while `config/paths.yml` uses the name for
`<deposit>/source_data`. Exporting it as the record root silently points those two scripts two
levels too high; they fail by finding nothing rather than by erroring. `DEPOSIT_ROOT` above avoids
the collision.

| Flag | Omitting it causes |
|---|---|
| `--containall` | The container loads the **host's** R packages and runs clean. Sites that set `mount home = yes` are not overridden by `--no-home` or `--cleanenv`. |
| `--pwd /work` | The container starts in an empty `$HOME`, and every repo-relative command fails with `does not exist` for a file that is plainly there. |
| `--nv` | **§5.7 only.** `torch.cuda.is_available()` returns False and §5 training runs on CPU **without an error**, by a different numerical path under mixed precision. Omit it for every other section. |

**`.libPaths()` cannot tell you whether `--containall` took.** Isolated and unisolated runs of the
same image report identical library paths. Assert isolation directly at the head of every stage:

```bash
[ -z "$(ls -A "$HOME" 2>/dev/null)" ] || { echo "NOT ISOLATED: \$HOME is not empty" >&2; exit 1; }
```

**Set `TMPDIR` and `SQLITE_TMPDIR`, and prefix both with `APPTAINERENV_`.** `--containall` gives
`/tmp` a 64 MB tmpfs and leaves `TMPDIR` unset. That is not enough for `productive_response.Rmd`'s
`AnnotationDbi::select()`, which fails with `database or disk is full` on a filesystem with
terabytes free. Point both at a bound directory with real space; setting the variables is
preferable to binding over `/tmp`, which leaves Apptainer's own tmpfs handling alone.

**`--containall` also discards the shell environment**, so a plain `export TMPDIR=…` arrives empty
inside the container and you hit the failure anyway. Apptainer passes a variable through only when
it is prefixed with `APPTAINERENV_`, and that applies to every variable you need inside — not only
the `NMD_*` ones REPRODUCTION.md names.

```bash
export SCRATCH=/path/with/real/space
mkdir -p "$SCRATCH/sqlite_tmp"
export APPTAINERENV_TMPDIR="$SCRATCH/sqlite_tmp"
export APPTAINERENV_SQLITE_TMPDIR="$SCRATCH/sqlite_tmp"

apptainer exec --containall \
  --bind "$REPO":/work --bind "$DEPOSIT_ROOT":"$DEPOSIT_ROOT" --bind "$SCRATCH":"$SCRATCH" \
  --pwd /work nmd_1.3.sif <command>
```

To confirm it took, run `echo $TMPDIR` inside the container: it should print the path, not an
empty line.

**Bind a directory, never a file.** Apptainer refuses a bind destination that does not exist,
where Docker creates one — so `--bind file:/tmp/x` works locally and fails in the image with a
mount error that reads like several unrelated defects.

## Versions

**R 4.5.2 (2025-10-31), Bioconductor 3.22**

**`renv.lock` is authoritative** and pins all 189. This table is a selection of the ones most
likely to move a number, read straight out of the lockfile — if the two ever disagree, the lockfile
is right.

| Package | Version | | Package | Version |
|---|---|---|---|---|
| data.table | 1.18.0 | | ashr | 2.2-63 |
| dplyr | 1.2.1 | | topGO | 2.62.0 |
| ggplot2 | 4.0.1 | | msigdbr | 26.1.0 |
| edgeR | 4.8.2 | | htmltools | 0.5.9 |
| tidyr | 1.3.1 | | readxl | 1.4.5 |
| tibble | 3.3.1 | | knitr | 1.51 |
| patchwork | 1.3.2 | | jsonlite | 2.0.0 |
| scales | 1.4.0 | | Biostrings | 2.78.0 |
| Isopair | 0.99.4 | | AnnotationDbi | 1.72.0 |
| reshape2 | 1.4.5 | | GenomicRanges | 1.62.1 |
| org.Hs.eg.db | 3.22.0 | | rtracklayer | 1.70.1 |
| mashr | 0.2.79 | | matrixStats | 1.5.0 |
| limma | 3.66.0 | | DT | 0.34.0 |

**`fgsea`, `pathview` and `tidyverse` were listed here with versions and are NOT in `renv.lock`.**
They are among the eight packages the container installs separately, so their versions are not
reproducible from the lockfile and must not be presented as if they were. They served code that
does not ship in this package, so they are simply removed rather than moved.

**Python 3.12.12**

| Package | Version | | Package | Version |
|---|---|---|---|---|
| matplotlib | 3.10.8 | | pillow (PIL) | 12.1.1 |
| numpy | 2.4.4 | | seaborn | 0.13.2 |
| pandas | 2.3.3 | | logomaker | 0.8.7 |
| scipy | 1.17.1 | | | |

**External tools**, as reported in Methods: nf-core/rnaseq v3.14.0, Nextflow 24.04.4, PacBio
Isocall, SQANTI3, TransDecoder2 (TD2), minimap2.

## Versions that change published numbers

Regenerating intermediates with different versions of these can change reported results.

| Package | Version | Affects |
|---|---|---|
| ORFik | 1.30.2 | `05s_orfik_scan.R` — the ORF scan behind the model's feature tables |
| msigdbr | 26.1.0 | GSEA gene-set membership (Tables 3–4) |
| mashr | 0.2.79 | Posterior means and lfsr (seeded, `set.seed(42)`) |
| edgeR / limma | 4.8.2 / 3.66.0 | Filtering, TMM normalization, differential testing |
| org.Hs.eg.db | 3.22.0 | §3.4 odds ratios |

**`org.Hs.eg.db` is data, not just a dependency.** `productive_response.Rmd` builds `THEME_SETS`
at run time from `org.Hs.egGO2ALLEGS` (`go_to_symbols`, line 158), so the published §3.4 odds
ratios are a function of this package's version: under 3.22.0, GO:0006260 resolves to 290 symbols
and GO:0006281 to 640. There is no file to deposit — the version pin is the provenance.

## Limitations

**ORFik.** 1.30.2 is the version in this environment. The originally deposited ORF scan predates
this record and its build version was not captured, so a regenerated scan is not guaranteed
byte-identical. If the scan drifts, Figure 5 and SF37–43 will not match the paper.

**Sections 1–3.** These results were produced in a separate R environment whose package versions
were not captured. The versions above are from the environment used for the remaining sections.
The versions that would matter most for §1–3 are `mashr`/`ashr` (posterior computation) and
`limma`/`edgeR` (the differential-expression fits).
