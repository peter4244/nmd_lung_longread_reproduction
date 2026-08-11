# Environment manifest — versions that produced the manuscript figures

Captured 2026-07-20 from Pete's laptop, the machine that rendered the tracked
figure PNGs. **This is a record, not a lockfile.** Full `renv::init()` is deferred to
the curated tree (see consolidation plan Phase 6.5 rationale below).

## Why a manifest now, `renv` later
Running `renv::init()` against the current repo would snapshot the dependencies of the
~200 files the prune is about to delete — locking in exactly the cruft being removed.
The irreplaceable information is *which versions produced the published figures*, so we
capture that now and generate the real lockfile from the curated tree at snapshot time.

## Corroboration
The tracked figure PNGs embed `Matplotlib version3.10.8`, matching the matplotlib below.
This is what makes the Phase 5.2 **PNG byte-compare** gate viable: the renderer on this
machine is the one that produced the committed figures.

> **CORRECTED 2026-08-10, and the original error is kept here because deleting it is how a
> wrong answer gets re-derived.** This section read *"embed `Matplotlib version3.11.0`"* and
> the table below pinned python 3.14.4, matplotlib 3.11.0, numpy 2.4.6, pandas 3.0.3, scipy
> 1.18.0 and pillow 12.2.0. Every one of those was a faithful capture — **of the wrong
> interpreter.** This laptop has two: Homebrew `python3.14`, which matches that set exactly,
> and miniforge `python3` at 3.12.12, which is what `python3 figure_*.py` actually runs. Of
> the 63 tracked PNGs, 36 are matplotlib's and **35 of those embed 3.10.8**; the single
> 3.11.0 was a container render that was reverted. So the corroboration sentence was
> asserting the opposite of what the artifacts say, and it was the stated justification for
> the pins — a restatement that had drifted from its source, which is the failure this
> project keeps paying for, sitting inside the document that records the environment.
>
> It surfaced only because the container rendered a panel that came out visibly different
> from the committed one, and the difference was measured rather than waved off as a
> platform quirk.

## R
```
R 4.5.2 (2025-10-31)   |   Bioconductor 3.22
```
| Package | Version | | Package | Version |
|---|---|---|---|---|
| data.table | 1.18.0 | | msigdbr | 26.1.0 |
| dplyr | 1.2.1 | | htmltools | 0.5.9 |
| ggplot2 | 4.0.1 | | readxl | 1.4.5 |
| edgeR | 4.8.2 | | knitr | 1.51 |
| tidyr | 1.3.1 | | jsonlite | 2.0.0 |
| tibble | 3.3.1 | | fgsea | 1.36.2 |
| patchwork | 1.3.2 | | Biostrings | 2.78.0 |
| scales | 1.4.0 | | AnnotationDbi | 1.72.0 |
| Isopair | 0.99.4 | | pathview | 1.50.0 |
| reshape2 | 1.4.5 | | GenomicRanges | 1.62.1 |
| org.Hs.eg.db | 3.22.0 | | rtracklayer | 1.70.1 |
| mashr | 0.2.79 | | matrixStats | 1.5.0 |
| limma | 3.66.0 | | tidyverse | 2.0.0 |
| hexbin | ⚠ NOT INSTALLED | | | |
| ashr | 2.2.63 | | DT | 0.34.0 |

**`topGO` 2.62.0 — installed locally as of 2026-07-28.** Used by
`analysis/upstream/productive_response.Rmd`. The earlier note here said it was NOT installed
and that §3 was parked; both were stale, and the note was found wrong while walking claim
3.4.5's path. Its two chunks are guarded by `eval = requireNamespace("topGO", ...)`, so on a
machine without it they silently do not run — which is why the absence went unnoticed.

**`org.Hs.eg.db` 3.22.0 is DATA, not just a dependency.** `productive_response.Rmd` builds
`THEME_SETS` at run time from `org.Hs.egGO2ALLEGS` (`go_to_symbols`, :158), so the published
odds ratios in §3.4 are a function of this package's version: under 3.22.0, GO:0006260 resolves
to 290 symbols and GO:0006281 to 640. There is no file to deposit and nothing for the provenance
graph to trace — the version pin above IS the provenance. Treat a change to it as a change to
the data (D25).

## Python
```
Python 3.12.12
```
| Package | Version |
|---|---|
| matplotlib | 3.10.8 |
| numpy | 2.4.4 |
| pandas | 2.3.3 |
| scipy | 1.17.1 |
| pillow (PIL) | 12.1.1 |
| seaborn | 0.13.2 |
| logomaker | 0.8.7 |

## External tools (pinned in Methods text; versions not yet verified against runs)
nf-core/rnaseq v3.14.0 · Nextflow 24.04.4 · PacBio Isocall · SQANTI3 · TransDecoder2 (TD2)
· minimap2. **Action (plan 6.5):** confirm each against the actual run logs, and cache the
`pathview` KEGG `hsa04141` template in-repo — it is currently fetched from the network at
runtime and will drift.

## The container, and the INVOCATION, which is part of the environment

The image pins R, the 187 packages in `renv.lock`, and both Python environments. **But an image
that is run wrongly is a different environment from the same image run correctly, so the
invocation is documented here rather than left in a job script.** A reader who obtains the image
and runs it the obvious way gets a different answer than we did — twice over, below.

**Where the image comes from: two ways.**

**Build it.** The recipe is `Dockerfile` in the repository root — it restores all 187 packages in
`renv.lock` and both conda environments, so it either succeeds or fails loudly. (The model
environment was gated behind `--build-arg WITH_MODEL=true`, which no documented build command
passed, and its else branch echoed and exited 0 — so until 2026-08-12 this path quietly produced
an image that could not run section 5. The default is now `true`.) This path works
today and is the one a reader can rely on.

**Or fetch it.** The Zenodo record carries `nmd_1.2.sif`, 8.0 GB — the image every result in the
paper came from, and the one the clean-room reruns used.

**Use `1.2`, not an older image.** `1.0` differs from `1.2` in two ways that reach the output: it
has **no Liberation Sans**, so figures fall through to DejaVu Sans, which is wider — the failure
that made SF42 pass the layout validator on the authoring laptop and fail in clean-room job
9049007 — and it predates the figures-environment repin, carrying matplotlib 3.11.0-era pins
rather than the 3.10.8 the tracked PNGs were actually rendered with.

```bash
docker build -t nmd:1.0 .
```

**Build and convert.** `docker save nmd:1.0 -o nmd.tar` then
`apptainer build nmd.sif docker-archive://nmd.tar`. Save a *locally built* image: a `docker save`
of a **pulled** multi-arch image produces a multi-manifest tarball Apptainer refuses, and adding
`--provenance` or SBOM attestations reintroduces that.

**`--pwd` and `--nv` are mandatory too, and neither names itself when omitted.** `--containall`
starts the container in an empty `$HOME` rather than the host working directory, and every command
in `REPRODUCTION.md` is repo-relative — so without `--pwd` the first step reports
`The file 'analysis/upstream/ct_de.Rmd' does not exist` while the file is plainly there. Without
`--nv` the container sees no GPU, `torch.cuda.is_available()` returns False, and the §5 training
step runs on CPU **without an error**, emitting a plausible checkpoint by a different numerical
path under mixed precision. The full form:

```bash
apptainer exec --containall --nv --bind "$REPO":/work --bind "$DEPOSIT":/deposit --pwd /work \
  nmd_1.2.sif <command>
```

**`--containall` is mandatory, not preference.** The site sets `mount home = yes` in
`/etc/apptainer/apptainer.conf`, which `--no-home --cleanenv` do **not** override, and Explorer
carries `~/R/x86_64-pc-linux-gnu-library`. Without `--containall` the container can load the
*host's* R packages, run clean, and report nothing.

**And `.libPaths()` cannot detect that it was omitted.** Measured 2026-08-05: isolated and
unisolated invocations of the same image both report
`/usr/local/lib/R/site-library`, `/usr/local/lib/R/library`. R skips the host library only because
that library is built for R 4.2 against the image's 4.5.2 — a version coincidence, not a barrier.
So **assert isolation at the head of every stage** rather than assuming it: check that `$HOME` is
*empty*, which is what actually proves `--containall` took, and that no authoring tree is
reachable. The script that does this in our working repository is not part of this package; the
check itself is two lines:

```bash
[ -z "$(ls -A "$HOME" 2>/dev/null)" ] || { echo "NOT ISOLATED: \$HOME is not empty" >&2; exit 1; }
ls /path/to/any/authoring/tree >/dev/null 2>&1 && { echo "NOT ISOLATED: host tree reachable" >&2; exit 1; }
```

**`--containall` also gives `/tmp` a 64 MB tmpfs, and leaves `TMPDIR` unset.** Measured: a 400 MiB
write into `/tmp` stops at exactly 64 MiB. That ceiling is enough for most of the chain and **not**
enough for `productive_response.Rmd`'s `AnnotationDbi::select()`, which fails with
`database or disk is full` while the cluster filesystem has hundreds of terabytes free. Point
`TMPDIR` **and** `SQLITE_TMPDIR` at a bound directory with real space:

```
--bind <scratch>:<scratch>   with  TMPDIR=<scratch>/sqlite_tmp_$SLURM_JOB_ID  SQLITE_TMPDIR=same
```

Setting the variables is preferred over binding a directory onto `/tmp`: it leaves Apptainer's own
tmpfs handling alone. Both work, and both keep the isolation assertions passing, because isolation
rides on `$HOME` being empty rather than on `/tmp`.

**Bind a DIRECTORY, never a file.** Apptainer refuses a bind destination that does not exist where
Docker creates one, so `--bind file:/tmp/x` works locally and fails in the image with a mount error
that reads like several unrelated defects.

**Verify the image itself** before trusting a run: check that R, both conda environments and
pandoc are present and at the expected versions. The script that does this in our working
repository is not part of this package —
the lockfile's version, all 187 packages at their pinned versions, `hexbin` and `mclust` (absent
from `renv.lock`), pandoc, and both Python environments.

## Not captured here
Yul's Channing environment (`/udd/reyle/Rlibs`), which is what actually produced the
§1–3 results (her scripts set `.libPaths()` to it). The versions above are Pete-side and
may differ. **Secondary ask** — worth requesting `sessionInfo()` from the manuscript run
when we next email her about the data bundle, but not a blocker on its own. The versions
that would actually matter are `mashr` / `ashr` (posterior computation) and `limma` /
`edgeR` (the DE fits), not the downstream annotation packages.

## Pinned package versions that affect published results

Regenerating intermediates with different versions of these can change published numbers, so
they are recorded explicitly:

| package | version | affects |
|---|---|---|
| ORFik | 1.30.2 | `05s_orfik_scan.R` — the ORF scan behind the model's feature tables |
| msigdbr | 26.1.0 | GSEA gene-set membership (Tables 3–4) |
| mashr | 0.2.79 | posterior means / lfsr (seeded, `set.seed(42)`) |
| edgeR / limma | 4.8.2 / 3.66.0 | filtering, TMM normalisation, differential testing |

⚠️ **ORFik caveat.** 1.30.2 is the version installed in the environment used for this work. The
originally deposited ORF scan predates this record and its exact build version was not captured,
so a regenerated scan is not guaranteed byte-identical to the original. A clean-room
reproduction is the check: if the scan drifts, Figure 5 and SF37–43 will not match the paper.

## Missing dependency — `hexbin`

`analysis/upstream/Figures/render_sr_lr_correlation.R` calls `stat_binhex()`, but `hexbin` is
not installed in this environment. The script still **exits 0** and writes its PNG — the hex
layer is simply absent, with only a warning. That is a silent-failure mode: a reproducer gets a
plausible-looking but empty panel and no error.

```r
install.packages("hexbin")
```

Found 2026-07-25 while repointing the SR↔LR scripts to the Zenodo deposit.
