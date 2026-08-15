# Build notes

For whoever rebuilds the container image. Reproducing the paper does not require this — fetch
`nmd_1.3.sif` from the Zenodo record instead. See [`ENVIRONMENT.md`](ENVIRONMENT.md).

These notes record why the `Dockerfile` is shaped the way it is, and the traps that cost real time
to find. Each is a measured finding, not a precaution.

## Build locally, convert on the cluster

Rootless image building is blocked on Explorer by site configuration: podman has no `/etc/subuid`
range, so even with the vfs driver on local disk it dies at `lchown /etc/shadow: invalid argument`;
apptainer gets further and fails with `exec /.singularity.d/libs/fakeroot failed`. Neither is
fixable without an administrator. Apptainer can *consume* an image perfectly well, so build here
and convert there:

```bash
docker build --platform linux/amd64 -t nmd:1.0 .
docker save nmd:1.0 -o nmd.tar
apptainer build nmd.sif docker-archive://nmd.tar     # on Explorer, no root needed
```

Three findings about that conversion, none visible from the Dockerfile itself:

1. **`--no-home --cleanenv` does not isolate `$HOME` on Explorer.** `/etc/apptainer/apptainer.conf`
   sets `mount home = yes` at site level, which those flags do not override. Explorer has a
   `~/R/x86_64-pc-linux-gnu-library` that would silently shadow the image's packages. Use
   `--containall` with explicit `--bind`; see `ENVIRONMENT.md`, which owns the run-time flags.
2. **`/tmp` on the login node does not persist** between an `scp` and a later `ssh`, so it is not a
   safe staging path for the tarball. Shared `/tmp` on Explorer has also collided with another
   user's file.
3. **`docker save` of a *pulled* multi-arch image yields a multi-manifest tarball Apptainer
   refuses** — `must specify a digest - layout contains multiple images`. A locally built image
   saves with one manifest and converts cleanly. Adding `--provenance` or SBOM attestations
   reintroduces extra manifests and breaks the conversion, and that failure surfaces much later
   looking unrelated.

## Platform

**`linux/amd64` explicitly.** The authoring laptop is arm64 and Explorer is x86_64. Left to
default, a build here produces an arm64 image that runs on neither.

**On Apple Silicon, use Rosetta, not qemu — qemu miscompiles.** With Colima on `vmType: vz` and
`rosetta: false`, amd64 emulation falls back to qemu binfmt and gcc died with
`Segmentation fault (core dumped)` compiling `globals.o` for vctrs. That is gcc itself crashing,
not a package failing to build, and any of the ~33 source builds can hit it — so the symptom moves
around and reads like a flaky package rather than a broken toolchain.

```bash
colima stop && colima start --vz-rosetta
```

Docker Desktop has the equivalent under Settings → General. A native x86_64 host needs none of it.

**Memory.** `Ncpus = 1` is what makes this build survivable on a small host: four concurrent C++
compiles exhausted an 8 GB VM. That is why the build is slow rather than parallel. On a machine
with real memory, raise it.

**The micromamba stage is platform-pinned for a reason.** An unpinned `COPY --from` resolves to the
*build host's* architecture, so on an ARM machine it copied an arm64 binary into an amd64 image and
the next step died with `micromamba: not found` — the file is present and is a perfectly good
22 MB ELF; it is the wrong ELF, and `not found` is what exec reports for that. Two earlier builds
succeeded with this defect fully present because they found an amd64 micromamba already in the
local image cache. **A successful build is not evidence the recipe is correct when the cache is
warm.**

## R packages, P3M, and a known unfixed defect

The first build compiled every package from source, because `renv::restore()` installs from the
repositories recorded in the **lockfile**, not from whatever the image has configured — and
`renv.lock`'s only repository is `cloud.r-project.org`, which serves source on Linux.

Two independent changes address that, and they fix different things:

1. **P3M serves prebuilt amd64 binaries for Ubuntu noble**, this image's platform. The `latest`
   endpoint was not usable — it had already moved past our stringi pin, which would send the
   worst offender straight back to source.
2. **`Ncpus = 1`** serialises whatever still compiles, so each fallback gets the whole VM rather
   than a quarter of it. This is the change that actually prevents the OOM; the binaries only make
   it fast enough to be worth doing.

**The snapshot date does not affect which versions are installed**, and that property is what makes
it safe. renv installs the version named in the lockfile: where the snapshot offers it, it is taken
as a binary; where it does not, renv falls back to source from the CRAN archive. The date is a
build-speed lever only. `2026-03-10` was chosen by measurement — it beat `2026-06-01` on exact
matches. No date serves all three heavy packages, so data.table remains a source build.

**P3M goes in front of the original repositories, never instead of them.** An earlier attempt set
`repos` to P3M alone, and a pin *newer* than the snapshot became unreachable: a dated snapshot has
no future packages and P3M's `Archive/` does not carry them. renv aborted the whole restore on one
package. `cloud.r-project.org` and `bioconductor.org` therefore stay in the list, last.

> ### Known defect, not yet fixed
>
> **The binary repository is largely not being used.** renv resolves a package against the repo
> whose *name* matches the `Repository` field recorded for that package in the lockfile. Most
> packages record `Repository: CRAN` — and the repos list names the P3M snapshot `P3M` while giving
> the name `CRAN` to `cloud.r-project.org`, which serves source on Linux. So renv goes to source
> for packages whose matching binaries sit unused in the snapshot.
>
> The projection that predicted otherwise was not wrong about **availability** — those versions are
> in the snapshot. It was wrong about **behaviour**: nothing checked that renv would look there.
>
> **The fix:** name the P3M snapshot `CRAN` so it wins the name match, and demote
> `cloud.r-project.org` to a differently-named entry (e.g. `CRANarchive`) that still serves as the
> archive fallback — renv searches all repos for an exact version when the named one lacks it. Same
> for Bioconductor: the P3M mirrors should carry the canonical `BioCsoft` / `BioCann` / `BioCexp`
> names.
>
> It was left alone because the restore layer was already cached at 6.68 GB and editing that `RUN`
> discards it — roughly 50 minutes of recompilation to save 40. **Apply it on the next rebuild that
> is needed for some other reason.**

## Packages installed outside renv.lock

Eight packages are installed by a separate `RUN` because they are loaded by shipped code with
`library()` and appear in no environment file. renv had no version to record, so the lockfile does
not mention them — and **the completeness guard only checks what `renv.lock` names**, so the
restore passed while this code still could not run.

| Package | Needed by | Note |
|---|---|---|
| `hexbin` | `correlation_analysis.Rmd`, two figure scripts | **Required.** Without it `geom_hex()` drops the layer, exits 0, and leaves a plausible empty panel |
| `mclust` | `productive_response.Rmd` | Genuinely optional — guarded by `requireNamespace`, so it silently skips when absent |
| `fgsea`, `pathview` | mashr summaries, pathway figures | Bioconductor |
| `tidyverse`, `xgboost`, `openxlsx` | DIE summaries, final report, table builder | |
| `magick` | `Figures/make_panels.R` | Needs `libmagick++-dev`; without it magick fails to *configure* and the R error names a header rather than a Debian package |

Their versions are **not reproducible from `renv.lock`** and must not be described as such:
the package was never present at publication, so what is installed here is the pinned snapshot's
current version.

## Why each guard exists

- **Restore completeness.** `renv::restore()` can report success while individual packages were
  skipped; a missing one would otherwise surface inside whichever `.Rmd` loads it first.
- **Python import check.** A solve can succeed, an install can succeed, and the result still be
  unimportable. One build pinned a `shap`/`numpy` combination where both resolved, both installed,
  `docker build` exited 0, and `import shap` raised
  `TypeError: Converting np.inexact or np.floating to a dtype not allowed`. The image built cleanly
  and could not run the two steps the model environment exists for.
- **Bare `Rscript` check, and it must stay after the conda environments are created.** One
  successful build shipped an image where bare `Rscript` resolved to conda's R, whose library held
  `pROC` and none of the analysis packages, while the system R held the analysis packages and not
  `pROC`. Neither could render the model report, and a reader got
  `there is no package called 'Isopair'` on a machine where it was plainly installed. The R restore
  guard ran *before* the conda environments existed, and therefore before `PATH` changed under it.
  The check uses bare `Rscript` deliberately: that is what `REPRODUCTION.md` tells readers to type.

## Fonts

`fonts-liberation` is installed by the apt layer, and it is the one thing the recipe originally did
not provide. The image the published figures were rendered under, `nmd_1.2.sif`, was **not** built
from this file — `apptainer inspect` reports it as `bootstrap: localimage, from: nmd_1.1.sif`, an
Apptainer layer over 1.1 whose whole delta is four Liberation Sans TTFs, with no Docker path to it.
The recipe installed font *libraries* and no font, so anyone who built rather than downloaded got a
container rendering every figure in whatever matplotlib fell back to.

Installing via apt rather than the original tarball is deliberate: the package puts the same four
files at the same path, and a version resolved from a pinned base image is more reproducible than a
tarball on a scratch filesystem. It may not be byte-identical to 1.2's, which is one more reason
figures are verified by comparing data exports rather than pixels.
