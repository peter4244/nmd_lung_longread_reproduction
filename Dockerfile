# The environment this paper's analyses run in — built from renv.lock and the two conda
# environment files, which are the recipes; this is the thing that proves they work.
#
# WHY THIS FILE EXISTS AT ALL. REPRODUCTION.md says a lockfile that has never been restored
# is an assertion. Building this image restores all 187 R packages and either succeeds or
# fails loudly, which is the whole point: it is the first guard in this repository that can
# actually fail on an environment claim. W292 stays open until it has been built and the
# chain has run inside it.
#
# WHY IT IS BUILT LOCALLY AND NOT ON THE CLUSTER. Rootless image building on Explorer is
# blocked by site configuration, measured 2026-08-03: podman has no /etc/subuid range, so
# even with the vfs driver on local disk it dies at `lchown /etc/shadow: invalid argument`;
# apptainer gets further and then fails with `exec /.singularity.d/libs/fakeroot failed`.
# Neither is fixable without an administrator. Apptainer can CONSUME an image perfectly well
# though, so building here and converting sidesteps both:
#
#     docker build --platform linux/amd64 -t nmd:1.0 .
#     docker save nmd:1.0 -o nmd.tar
#     apptainer build nmd.sif docker-archive://nmd.tar     # on Explorer, no root needed
#
# THE CONVERSION IS MEASURED, NOT ASSUMED (2026-08-05, on Explorer). Unprivileged
# `apptainer build` from a docker archive works there. Three site-level findings that cost real
# time to discover and are invisible from the Dockerfile alone:
#
#   1. `--no-home --cleanenv` DOES NOT ISOLATE $HOME on Explorer. /etc/apptainer/apptainer.conf
#      sets `mount home = yes` at SITE level, which those flags do not override. This matters
#      because Explorer has a ~/R/x86_64-pc-linux-gnu-library that would silently SHADOW the
#      image's packages -- the pinned environment defeated by a stray user library, with no
#      error. A flag that looks like it prevents this and does not is worse than no flag.
#      Use `--containall` plus explicit `--bind` for what the pipeline actually needs; both the
#      isolation and the required reads/writes were verified.
#   2. /tmp on the LOGIN NODE does not persist between an scp and a later ssh. It is not a safe
#      staging path for the tarball. (Separately: never use shared /tmp on Explorer at all --
#      a redirect there once collided with another user's file.)
#   3. `docker save` of a PULLED multi-arch image yields a multi-manifest tarball that Apptainer
#      refuses: "must specify a digest - layout contains multiple images". A LOCALLY BUILT image
#      saves with exactly one manifest and converts cleanly, so nmd:1.0 is fine -- but adding
#      `--provenance` or SBOM attestations to the build reintroduces extra manifests and breaks
#      the conversion. That failure would surface much later and look unrelated to this file.
#
# WHY linux/amd64 EXPLICITLY. The authoring laptop is arm64 and Explorer is x86_64. Left to
# default, a build here produces an arm64 image that cannot run there and that most readers
# cannot run either. Building amd64 on an arm64 host uses emulation, which is slow for
# COMPILATION but cheap for downloads — which is why the binary package repository below
# matters more than it looks.
#
# IF YOU BUILD THIS ON APPLE SILICON, USE ROSETTA, NOT QEMU. This is not a performance
# preference; qemu MISCOMPILES. Measured here 2026-08-05: with Colima on `vmType: vz` and
# `rosetta: false`, amd64 emulation falls back to qemu binfmt and gcc died with
# `Segmentation fault (core dumped)` compiling globals.o for vctrs 0.7.2 — gcc itself crashing,
# not the package failing to build. Any of the ~33 source builds can hit it, so the symptom
# moves around and reads like a flaky package rather than a broken toolchain.
#
#     colima stop && colima start --vz-rosetta        # then rebuild
#
# Docker Desktop has the equivalent switch under Settings > General ("Use Rosetta for
# x86_64/amd64 emulation"). A native x86_64 host needs none of this.
#
# THE HOST ALSO NEEDS MEMORY, AND THE AUTHORING LAPTOP BARELY HAD IT. 8 GB total RAM, all 8
# allocated to the VM, is where the FIRST failure came from: four concurrent C++ compiles
# exhausted it. Ncpus = 1 below is what makes that survivable, and it is why this build is
# slow rather than parallel. On a machine with real memory, raise it.
#
# WHAT IS DELIBERATELY NOT IN THE IMAGE: the analysis code and the data deposit. The image is
# the environment and nothing else, so re-running does not mean rebuilding, and so the image
# cannot silently carry a stale copy of the code it is supposed to be neutral about. Both
# bind-mount at run time:
#
#     docker run --rm -v "$PWD":/work -v /path/to/nmd_deposit_2026:/deposit nmd:1.0 \
#            python3 tools/run_chain.py --room /work/room --run-id docker-1 --select ""

ARG TARGETPLATFORM=linux/amd64

# R 4.5.2 EXACTLY, because renv.lock pins it and renv warns on a mismatch. rocker/r-ver is
# preferred over the Bioconductor image for precisely that reason: the Bioconductor image is
# keyed to a release (3.22) and tracks whatever R patch that release carries, which is not a
# guarantee of 4.5.2.
#
# THIS COMMENT USED TO CLAIM rocker's Posit binary repository would make most packages install
# as prebuilt binaries. THAT IS FALSE HERE, and believing it cost a build. renv::restore()
# installs from the repositories recorded in the LOCKFILE, not from whatever the image has
# configured, and renv.lock's only repository is https://cloud.r-project.org -- which serves
# SOURCE on Linux. So every one of the 182 packages compiled, four at a time, and stringi,
# data.table and RcppParallel (the three heaviest C++ compiles in the set) all died at once on
# a 4-CPU / 7.7 GiB VM emulating amd64. See the P3M block below for the actual fix, and
# /private/tmp/nmd_docker_build_attempt1_OOM.log for the failure.
FROM --platform=linux/amd64 rocker/r-ver:4.5.2

ARG DEBIAN_FRONTEND=noninteractive
ARG PANDOC_VERSION=3.8.3
# Set to true to add the model environment (torch). Off by default: it pulls a CUDA build of
# torch, roughly 3 GB, which is useless without an NVIDIA GPU. Sections 1-4 do not touch it.
# DEFAULT true SINCE 2026-08-12. It was false, and no documented build command overrode it, so
# every build the documentation describes produced an image WITHOUT the model environment -- while
# README.md and ENVIRONMENT.md both promised both environments and claimed the build "fails loudly
# if it cannot". The else branch below echoes and exits 0, so it failed silently and section 5 then
# had no interpreter. Build with --build-arg WITH_MODEL=false to opt out deliberately.
ARG WITH_MODEL=true

# System libraries the Bioconductor packages link against. These are the ones that actually
# bite: rtracklayer and Biostrings need libxml2/zlib/bzip2/lzma/curl; the graphics stack
# (ragg, Cairo) needs freetype/harfbuzz/fribidi/png/jpeg/tiff; topGO pulls graph/Rgraphviz;
# glpk is igraph's. Missing any one of them surfaces as a COMPILER error partway through
# renv::restore(), not as a missing-package error, which is why they are installed up front
# rather than discovered one failed build at a time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran \
        libxml2-dev libcurl4-openssl-dev libssl-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        libpng-dev libjpeg-dev libtiff-dev \
        libcairo2-dev libfontconfig1-dev libfreetype6-dev \
        libharfbuzz-dev libfribidi-dev \
        libglpk-dev libgsl-dev libgit2-dev libxt-dev \
        git ca-certificates curl bzip2 \
        fonts-liberation \
    && rm -rf /var/lib/apt/lists/*

# fonts-liberation IS ABOVE, AND IT IS THE ONE THING THE RECIPE DID NOT BUILD.
#
# The image the clean room actually ran -- nmd_1.2.sif, which produced job 9057341 -- was NOT built
# from this file. `apptainer inspect` reports it as `bootstrap: localimage, from: nmd_1.1.sif`: an
# Apptainer-side layer over 1.1, with no Docker path to it at all. Its whole 860 KB delta is four
# Liberation Sans TTFs, extracted from a tar on /scratch by a def file that lived only there. The
# recipe installed font LIBRARIES (libfontconfig1-dev, libfreetype6-dev) and NO FONT, so a reader
# following ENVIRONMENT.md's "you build it ... there is no image to download" got a container that
# renders every figure in whatever matplotlib falls back to.
#
# Measured 2026-08-11 on Explorer: nmd_1.1.sif carries 0 Liberation fonts, nmd_1.2.sif carries 4,
# matplotlib inside 1.2 resolves "Liberation Sans", and Arial is absent from both. The def file is
# preserved at tools/cleanroom/nmd_1.2_fonts.def because /scratch is not backed up and it was the
# only record of how the image we ran came to exist.
#
# APT RATHER THAN THE TAR, deliberately: the package puts the same four files at the same path
# (/usr/share/fonts/truetype/liberation/LiberationSans-{Regular,Bold,Italic,BoldItalic}.ttf), and a
# version resolved from a pinned base image is more reproducible than a tarball on a scratch
# filesystem. It may not be byte-identical to 1.2's tar -- which is one more reason the pixel gate
# stays off under D120.

# pandoc at the PINNED version, not the distribution's. Ubuntu ships a much older pandoc;
# rmarkdown would accept it (it needs >= 1.12.3) and the reports would render, but they would
# render under a different pandoc than the one that produced ours, which is exactly the kind
# of undeclared difference this repository keeps getting caught by. pandoc is also invisible
# to renv, so nothing else in the reproducibility apparatus would notice.
RUN curl -fsSL -o /tmp/pandoc.deb \
        "https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-amd64.deb" \
    && dpkg -i /tmp/pandoc.deb && rm /tmp/pandoc.deb \
    && pandoc --version | head -1

# --- R packages -------------------------------------------------------------------------
# Copied alone, before anything else that changes often, so edits elsewhere do not invalidate
# the most expensive layer in the build.
COPY renv.lock /tmp/renv.lock

# renv::restore() into the site library rather than a project library: there is no renv
# project inside the image, and the code that runs here is bind-mounted from outside, so a
# project-scoped library would not be on its path.
#
# Isopair comes from GitHub (peter4244/Isopair at ea33f74) and is fetched by this step, which
# is why git and ca-certificates are installed above.
#
# GITHUB_PAT IS NOT DEFAULTED TO "", AND THAT COMMENT USED TO SAY THE OPPOSITE. An adversarial
# review reproduced the failure end to end: ARG values are exported into the environment of
# later RUN steps, so `ARG GITHUB_PAT=""` sets the variable to an EMPTY STRING rather than
# leaving it unset. renv tests `Sys.getenv(envvar, unset = NA)` and takes any non-NA value, so
# it emits the header `Authorization: token ` with nothing after it -- and GitHub answers 401
# Bad credentials on a PUBLIC repo that downloads fine unauthenticated. Isopair has a RemoteSha
# but no RemoteUrl, so renv registers no git fallback for it and there is no recovery path:
# `failed to install Isopair`, and the build dies.
#
# So the empty default was not a rate-limit convenience, it was the thing that broke the build.
# Unset it in the same shell unless a real token was passed.
#
# A real token would also be baked into image metadata and would ship inside the .tar and the
# .sif. If one is ever genuinely needed, use a BuildKit secret mount rather than --build-arg.
# PREBUILT BINARIES FROM A DATE-PINNED POSIT SNAPSHOT, AND WHY THE DATE IS NOT ARBITRARY.
#
# The first build compiled all 182 packages from source because renv uses the lockfile's
# repositories (CRAN, source-only on Linux) and ignores the image's. Four concurrent C++
# compiles under qemu on a 7.7 GiB VM exhausted memory and killed stringi, data.table and
# RcppParallel. Two independent changes below, because they fix two different things:
#
#   1. P3M serves prebuilt amd64 binaries for Ubuntu noble, which is exactly this image's
#      platform. Measured against renv.lock at this snapshot: 153 of 187 packages match the
#      pinned version exactly and install WITHOUT COMPILING; 30 differ in version; 3 are absent.
#      The `latest` endpoint was NOT usable -- it had already moved to stringi 1.8.9 against our
#      1.8.7 pin, which would have sent the worst offender straight back to source.
#
#   2. Ncpus = 1 serialises whatever still has to compile, so each of the ~33 fallbacks gets the
#      whole 7.7 GiB rather than a quarter of it. This is the change that actually prevents the
#      OOM; the binaries just make it fast enough to be worth doing.
#
# THE SNAPSHOT DATE DOES NOT AFFECT WHICH VERSIONS ARE INSTALLED, and that property is the
# reason this is safe. renv installs the version named in the lockfile: where the snapshot
# offers that version it is taken as a binary, and where it does not renv falls back to source
# from the CRAN archive. The date is a build-speed lever ONLY. A reader who changes it gets a
# faster or slower build, never a different library -- and the verification step that follows
# would catch it if that were ever untrue.
#
# 2026-03-10 was chosen by measurement, not by intuition: it beat 2026-06-01 by 153 exact
# matches to 130. No date serves all three heavy packages (data.table 1.18.0 needs January,
# RcppParallel 5.1.11-2 needs March), so data.table remains a source build and is fine alone.
#
# P3M GOES IN FRONT OF THE ORIGINAL REPOSITORIES, NEVER INSTEAD OF THEM. Attempt 2 failed here:
# it set `repos` to P3M alone, and `msigdbr 26.1.0` -- a pin NEWER than the snapshot -- became
# unreachable, because a dated snapshot has no future packages and P3M's Archive/ does not carry
# them either. renv tried four P3M archive paths, found nothing, and aborted the whole restore on
# one package. The snapshot is an ACCELERATOR layered over the real repositories; the moment it
# is treated as a replacement, every pin outside its date window has nowhere to come from.
# cloud.r-project.org and bioconductor.org therefore stay in the list, last, as the fallback that
# actually holds the archive. Order is what makes this both fast and complete: renv takes the
# first repository that has the exact pinned version, so binaries win where they exist and source
# is always reachable where they do not.
# THE BINARY REPOSITORY BELOW IS LARGELY NOT BEING USED, AND THE REASON IS A DEFECT IN THIS
# FILE. Measured on the 2026-08-05 build: 32 of 182 packages installed as binaries, not the
# 153/187 the snapshot predicted. renv resolves a package against the repo whose NAME matches
# the `Repository` field recorded FOR THAT PACKAGE in the lockfile. 147 packages record
# `Repository: CRAN` -- and the repos list below names the P3M snapshot `P3M` while giving the
# name `CRAN` to cloud.r-project.org, which serves SOURCE on Linux. So renv went to source for
# 142 packages with the matching binaries sitting unused in the snapshot.
#
# The projection was not wrong about AVAILABILITY -- those exact versions are in the snapshot,
# which was verified against renv.lock before the build. It was wrong about BEHAVIOUR: nobody
# checked that renv would look there. Availability is not behaviour, and only the first was
# tested.
#
# THE FIX, DELIBERATELY NOT APPLIED MID-BUILD: name the P3M snapshot `CRAN` so it wins the name
# match, and demote cloud.r-project.org to a differently-named entry (e.g. `CRANarchive`) that
# still serves as the archive fallback -- renv searches ALL repos for an exact version when the
# named one lacks it, which is what attempt 2 demonstrated when it walked four archive URLs
# before failing on msigdbr. Same for Bioconductor: the P3M mirrors should carry the canonical
# `BioCsoft` / `BioCann` / `BioCexp` names and bioconductor.org should take the fallback names.
#
# It was left alone on 2026-08-05 because the restore layer was already cached at 6.68 GB, and
# editing this RUN discards it: ~50 minutes of recompilation spent to save ~40. Apply it on the
# next rebuild that is needed for some other reason.
ARG P3M_SNAPSHOT=2026-03-10
ARG BIOC_VERSION=3.22

ARG GITHUB_PAT=""
RUN [ -n "$GITHUB_PAT" ] || unset GITHUB_PAT; \
    P3M="https://p3m.dev/cran/__linux__/noble/${P3M_SNAPSHOT}"; \
    BIO="https://p3m.dev/bioconductor/__linux__/noble/packages/${BIOC_VERSION}"; \
    R -q -e "install.packages('renv', repos = '${P3M}')" \
 && R -q -e "options(Ncpus = 1, \
                     repos = c(P3M       = '${P3M}', \
                               P3MBioc   = paste0('${BIO}', '/bioc'), \
                               P3MBiocAnn= paste0('${BIO}', '/data/annotation'), \
                               P3MBiocExp= paste0('${BIO}', '/data/experiment'), \
                               CRAN      = 'https://cloud.r-project.org', \
                               BioCsoft  = 'https://bioconductor.org/packages/${BIOC_VERSION}/bioc', \
                               BioCann   = 'https://bioconductor.org/packages/${BIOC_VERSION}/data/annotation', \
                               BioCexp   = 'https://bioconductor.org/packages/${BIOC_VERSION}/data/experiment')); \
             renv::restore(lockfile = '/tmp/renv.lock', library = .libPaths()[1], \
                           prompt = FALSE, repos = getOption('repos'))"

# Fail the BUILD, not some analysis step three hours into a run, if the restore was partial.
# renv::restore() can report success while individual packages were skipped, and a missing
# package would otherwise surface as an error inside whichever .Rmd happens to load it first.
RUN R -q -e 'pkgs <- names(renv::lockfile_read("/tmp/renv.lock")$Packages); \
             miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
             if (length(miss)) stop("restore incomplete: ", paste(miss, collapse = ", ")); \
             cat("all", length(pkgs), "packages load\n")'

# hexbin AND mclust, INSTALLED HERE AND ABSENT FROM renv.lock (W293, revised 2026-08-05).
#
# Both are referenced by keep-list code and installed on no machine we can inspect, so renv had
# no version to record and the lockfile does not mention them. An earlier version of this file
# left them out on the grounds that installing would "pin whatever version today happens to
# supply and assert a provenance we do not have."
#
# THAT REASONING WAS IMPORTED FROM A CASE IT DOES NOT FIT. It holds when a known version is
# being replaced by a guess. Here NO version was ever installed, so there is no original
# provenance to protect -- the choice is not true-version against false-version, it is some
# version against a silently broken result. Specifically:
#
#   hexbin   REQUIRED, not optional. geom_hex() in analysis/upstream/correlation_analysis.Rmd
#            (lines 324, 425) and two figure scripts. ENVIRONMENT.md has recorded it missing
#            since 2026-07-20 together with its failure mode: the layer drops, exit 0,
#            plausible empty panel. So every run to date has rendered those panels degraded
#            and nothing said so.
#   mclust   genuinely optional. productive_response.Rmd:1027 guards it with
#            eval = requireNamespace(...) and labels the chunk exploratory, so absent it
#            silently skips by design. Installed anyway so the chunk actually runs rather
#            than being permanently dark.
#
# Their versions are therefore NOT reproducible from renv.lock and must not be described as
# such. Recorded as: version unknown at publication because the package was never present;
# current CRAN installed here.
RUN R -q -e 'install.packages(c("hexbin", "mclust"), repos = "https://cloud.r-project.org"); \
             for (p in c("hexbin", "mclust")) if (!requireNamespace(p, quietly = TRUE)) \
               stop("failed to install ", p); \
             cat("hexbin", as.character(packageVersion("hexbin")), \
                 " mclust", as.character(packageVersion("mclust")), "\n")'

# --- Python -----------------------------------------------------------------------------
# micromamba as a single static binary from its own image — no installer script, no download
# URL to rot, and it resolves the two environment files exactly as conda would.
COPY --from=mambaorg/micromamba:latest /bin/micromamba /usr/local/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Both environments, kept separate for the reason environment-figures.yml gives: they ran
# under different Pythons, and merging them would mean choosing one matplotlib when
# matplotlib 3.11.0 is the pin the figure byte-comparison depends on.
COPY environment-figures.yml /tmp/environment-figures.yml
COPY environment-model.yml   /tmp/environment-model.yml

RUN micromamba create -y -f /tmp/environment-figures.yml && micromamba clean -a -y

RUN if [ "$WITH_MODEL" = "true" ]; then \
        micromamba create -y -f /tmp/environment-model.yml && micromamba clean -a -y; \
    else \
        echo "WARNING: model environment SKIPPED via --build-arg WITH_MODEL=false. Section 5 cannot run in this image."; \
    fi

# The chain rewrites the leading `python3`/`Rscript` of each step command through these, which
# is what lets one image serve steps that expect different interpreters (tools/run_chain.py,
# resolve_interpreters). Pointing PYTHON at the figures environment is the deliberate default:
# it is the one every non-model step needs.
ENV PYTHON=/opt/conda/envs/nmd_figures/bin/python \
    RSCRIPT=/usr/local/bin/Rscript \
    PATH=/opt/conda/envs/nmd_figures/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
