# The environment this paper's analyses run in — built from renv.lock and the two conda
# environment files, which are the recipes; this is the thing that proves they work.
#
# Building this image restores all 189 R packages in renv.lock and both conda environments, and
# either succeeds or fails loudly — a lockfile that has never been restored is an assertion.
#
# You do not need to build it to reproduce the paper. Fetch nmd_1.3.sif from the Zenodo record
# instead; see ENVIRONMENT.md.
#
#     docker build --platform linux/amd64 -t nmd:1.0 .
#     docker save nmd:1.0 -o nmd.tar
#     apptainer build nmd.sif docker-archive://nmd.tar     # on a cluster, no root needed
#
# THREE THINGS THAT WILL BITE YOU, each measured. BUILD_NOTES.md has the detail.
#
#   1. Build LOCALLY and convert on the cluster — rootless image building is blocked by site
#      configuration on Explorer, and is not fixable without an administrator.
#   2. `docker save` of a PULLED multi-arch image yields a multi-manifest tarball Apptainer
#      refuses. Save a LOCALLY BUILT image; `--provenance` and SBOM attestations reintroduce
#      the extra manifests.
#   3. ON APPLE SILICON, USE ROSETTA, NOT QEMU — qemu MISCOMPILES. gcc itself segfaults partway
#      through a source build, which reads like a flaky package rather than a broken toolchain.
#      `colima stop && colima start --vz-rosetta`, or Docker Desktop Settings > General.
#
# --platform linux/amd64 is explicit because an arm64 build runs on neither the cluster nor most
# readers' machines. Ncpus = 1 below is what keeps the build inside a small host's memory; raise
# it on a larger machine.
#
# The analysis code and the data deposit are deliberately NOT in the image — it is the
# environment and nothing else, so it cannot carry a stale copy of the code. Both bind-mount at
# run time:
#
#     docker run --rm -v "$PWD":/work -v /path/to/nmd_deposit_2026:/deposit nmd:1.0 <command>

ARG TARGETPLATFORM=linux/amd64

# MICROMAMBA AS A NAMED, PLATFORM-PINNED STAGE. Both pins are load-bearing and neither was here
# before 2026-08-12.
#
# --platform: this was `COPY --from=mambaorg/micromamba:latest` inline, with no platform. Every
# other image here is pinned to linux/amd64, but an unpinned COPY --from resolves to the BUILD
# HOST's architecture -- so on an ARM machine it copied an arm64 binary into an amd64 image, and
# the next step died with `micromamba: not found`. The file is present and is a perfectly good
# 22 MB ELF; it is the wrong ELF, and "not found" is what exec reports for that.
#
# THE PART WORTH RECORDING: this recipe built successfully on an ARM laptop on 2026-08-06 and
# again on 2026-08-10 with this defect fully present. Those builds found an amd64 micromamba
# already in the local image cache and copied that. Wiping the VM forced a fresh pull, the pull
# resolved to arm64, and the latent defect surfaced. So two successful builds were never evidence
# the recipe was correct -- they were evidence the cache was warm.
#
# :latest -> 2.9.0: an unpinned tag in a recipe whose entire purpose is reproducing an
# environment. The version that this pin was tested against is the one recorded here.
FROM --platform=linux/amd64 mambaorg/micromamba:2.9.0 AS micromamba

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
        libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*
# libmagick++-dev is for the R `magick` package, added 2026-08-12 with the six packages installed
# further down. Without it magick fails to CONFIGURE, and the R error names a header rather than a
# Debian package, which is a long way from "apt-get install libmagick++-dev".

# fonts-liberation IS LOAD-BEARING, and the recipe originally omitted it. It installed font
# LIBRARIES (libfontconfig1-dev, libfreetype6-dev) and NO FONT, so anyone who built this image
# rather than downloading one got a container rendering every figure in whatever matplotlib fell
# back to. Arial is licensed and cannot be redistributed; Liberation Sans is the metric-compatible
# substitute. Installed via apt so it resolves from the pinned base image. See BUILD_NOTES.md.

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
# PREBUILT BINARIES FROM A DATE-PINNED POSIT SNAPSHOT. Without this every package compiles from
# source: renv::restore() installs from the repositories recorded in the LOCKFILE, not from
# whatever the image has configured, and renv.lock's only repository serves SOURCE on Linux.
#
# Two independent changes, fixing two different things. P3M serves prebuilt amd64 binaries for
# Ubuntu noble, this image's platform. Ncpus = 1 serialises whatever still has to compile, so
# each fallback gets the whole VM rather than a quarter of it -- that is what prevents the OOM;
# the binaries only make it fast enough to be worth doing.
#
# THE SNAPSHOT DATE IS A BUILD-SPEED LEVER ONLY, which is what makes it safe. renv installs the
# version named in the lockfile: where the snapshot offers it, it is taken as a binary; where it
# does not, renv falls back to source from the CRAN archive. Changing the date gives a faster or
# slower build, never a different library.
#
# P3M GOES IN FRONT OF THE ORIGINAL REPOSITORIES, NEVER INSTEAD OF THEM. A dated snapshot holds
# nothing newer than its date, so treating it as a replacement leaves every pin outside that
# window with nowhere to come from and aborts the whole restore on one package.
# cloud.r-project.org and bioconductor.org stay in the list, last, as the archive fallback.
# !! KNOWN DEFECT, NOT YET FIXED: THE BINARY REPOSITORY BELOW IS LARGELY UNUSED. renv resolves a
# package against the repo whose NAME matches the `Repository` field recorded for that package in
# the lockfile. Most packages record `Repository: CRAN` -- and this list names the P3M snapshot
# `P3M` while giving the name `CRAN` to cloud.r-project.org, which serves SOURCE on Linux. So renv
# goes to source for packages whose matching binaries sit unused in the snapshot. It costs build
# time only; the installed versions are correct either way.
#
# The fix, and why it was deferred rather than applied mid-build, are in BUILD_NOTES.md. Apply it
# on the next rebuild needed for some other reason.
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

# EIGHT PACKAGES INSTALLED HERE AND ABSENT FROM renv.lock. They are loaded by shipped code with
# library() and appear in no environment file, so renv had no version to record. THE COMPLETENESS
# GUARD ABOVE CHECKS ONLY WHAT renv.lock NAMES, so the restore passes while this code cannot run.
#
#   hexbin     REQUIRED. correlation_analysis.Rmd and two figure scripts. Without it geom_hex()
#              drops the layer, exits 0, and leaves a plausible empty panel
#   mclust     optional -- productive_response.Rmd guards it with requireNamespace
#
# Their versions are NOT reproducible from renv.lock and must not be described as such: the
# package was never present at publication, so this installs the pinned snapshot's version.
#
#   fgsea      Gene-Level_DGE_Summary_mashR.Rmd, interpret_isoform_patterns_mashr   [Bioconductor]
#   pathview   Figures/make_pathway_allct.R                                          [Bioconductor]
#   tidyverse  Isoform-Level_DIE_Summary_p1.Rmd, die_mashr_enrichment_part2
#   xgboost    isopair_wrapper/05_final_report_mashr.Rmd
#   openxlsx   build_manuscript_tables.R
#   magick     Figures/make_panels.R  -- needs libmagick++-dev, added to the apt layer above
#
# P3M repos, not cloud.r-project.org, because two are Bioconductor and the CRAN mirror cannot
# resolve them. Same pinned snapshot the restore uses, so these land at the same vintage.
ARG P3M_SNAPSHOT_EXTRA=2026-03-10
RUN P3M="https://p3m.dev/cran/__linux__/noble/${P3M_SNAPSHOT_EXTRA}"; \
    BIO="https://p3m.dev/bioconductor/__linux__/noble/packages/${BIOC_VERSION}"; \
    R -q -e "pkgs <- c('hexbin','mclust','fgsea','pathview','tidyverse','xgboost','openxlsx','magick'); \
             install.packages(pkgs, repos = c(P3M = '${P3M}', \
                                              P3MBioc = paste0('${BIO}', '/bioc'), \
                                              P3MBiocAnn = paste0('${BIO}', '/data/annotation'), \
                                              CRAN = 'https://cloud.r-project.org', \
                                              BioCsoft = 'https://bioconductor.org/packages/${BIOC_VERSION}/bioc')); \
             miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
             if (length(miss)) stop('failed to install: ', paste(miss, collapse = ', ')); \
             for (p in pkgs) cat(p, as.character(packageVersion(p)), '\n')"

# --- Python -----------------------------------------------------------------------------
# micromamba as a single static binary from its own image — no installer script, no download
# URL to rot, and it resolves the two environment files exactly as conda would. The stage is
# declared at the top of this file with an explicit --platform; see the comment there for why
# copying from an unpinned image reference is not the same thing.
COPY --from=micromamba /bin/micromamba /usr/local/bin/micromamba

# Fail HERE if the binary cannot execute, rather than three steps later with `micromamba: not
# found`, which reads like a PATH problem and is not one.
RUN micromamba --version
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Both environments, kept separate for the reason environment-figures.yml gives: they ran
# under different Pythons, and merging them would mean choosing one matplotlib when the figure
# comparison depends on the pinned 3.10.8 — see the version note in environment-figures.yml.
COPY environment-figures.yml /tmp/environment-figures.yml
COPY environment-model.yml   /tmp/environment-model.yml

RUN micromamba create -y -f /tmp/environment-figures.yml && micromamba clean -a -y

RUN if [ "$WITH_MODEL" = "true" ]; then \
        micromamba create -y -f /tmp/environment-model.yml && micromamba clean -a -y; \
    else \
        echo "WARNING: model environment SKIPPED via --build-arg WITH_MODEL=false. Section 5 cannot run in this image."; \
    fi

# IMPORT WHAT WAS INSTALLED. A solve can succeed, an install can succeed, and the result can still
# be unimportable -- so an install check is not an environment check, exactly as renv::restore()'s
# exit code is not a check that 189 packages load. This exists because the file once
# pinned shap 0.46 against numpy 2.4: both resolved, both installed, `docker build` exited 0, and
# `import shap` raised `TypeError: Converting np.inexact or np.floating to a dtype not allowed`.
# The image built cleanly and could not run the DeepSHAP or KernelSHAP steps -- the two steps the
# model environment exists for.
#
# Import the packages the SHIPPED CODE actually imports, not a token one. Anything cheap to add
# here is cheaper than a reader discovering it four steps into a cluster run.
RUN /opt/conda/envs/nmd_figures/bin/python -c \
      "import matplotlib, numpy, pandas, scipy, seaborn, logomaker; print('figures env: imports OK')" \
 && if [ "$WITH_MODEL" = "true" ]; then \
        /opt/conda/envs/nmd_model/bin/python -c \
          "import torch, shap, h5py, numpy, pandas, yaml, tqdm; print('model env: imports OK')"; \
    fi

# AND THE SAME CHECK FOR R, RUN THROUGH BARE `Rscript` ON PURPOSE. Every command in
# REPRODUCTION.md is written as bare `Rscript`, so that is what has to work -- not the absolute
# path we happen to know. This is deliberately not `/usr/local/bin/Rscript`: resolving it the way
# the documentation does is the entire point of the check.
#
# It exists because on 2026-08-13 a SUCCESSFUL build shipped an image where bare `Rscript` resolved
# to conda's R, whose only library holds pROC and none of the 189 packages, while the system R held
# the 189 and not pROC. Neither could render the section 5 report, and a reader following the page
# got "there is no package called 'Isopair'" on a machine where it was plainly installed. Nothing
# in the build noticed, because the R restore guard ran BEFORE the conda environments existed and
# therefore before PATH changed under it.
#
# Ordering note: this must stay AFTER the conda environments are created, or it tests a PATH that
# the finished image will not have.
RUN Rscript -e 'need <- c("Isopair","ORFik","edgeR","mashr","pROC","ggseqlogo","rmarkdown"); \
                miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]; \
                if (length(miss)) stop("bare Rscript resolves to an R that is missing: ", \
                                       paste(miss, collapse = ", "), \
                                       " -- .libPaths() = ", paste(.libPaths(), collapse = " ; "), \
                                       " -- which Rscript = ", Sys.which("Rscript")); \
                cat("bare Rscript: resolves to the analysis library, all", length(need), "present\n")'

# The chain rewrites the leading `python3`/`Rscript` of each step command through these, which
# is what lets one image serve steps that expect different interpreters (tools/run_chain.py,
# resolve_interpreters). Pointing PYTHON at the figures environment is the deliberate default:
# it is the one every non-model step needs.
ENV PYTHON=/opt/conda/envs/nmd_figures/bin/python \
    RSCRIPT=/usr/local/bin/Rscript \
    PATH=/opt/conda/envs/nmd_figures/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
