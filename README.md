# Nonsense-mediated decay in human lung cells — analysis code

Code for the analyses and figures in *Long Read RNA Sequencing in Primary Lung Cell Types Reveals
Principles of Nonsense-Mediated Decay* (Leshem et al., 2026).

**The study in one paragraph.** Nonsense-mediated decay (NMD) is the cell's system for
destroying faulty messenger RNAs. We inhibited it in cultured human lung cells with an SMG1
inhibitor (`Smg1i` in the data files, against a `DMSO` vehicle control), sequenced the RNA with
long reads so that whole transcripts could be read end to end, and compared what appeared when
NMD was inhibited against what appeared in untreated cells. That tells us which RNA molecules
the cell had been quietly destroying. We did this in four lung cell types — alveolar type 2
cells (`AT2`), airway epithelial cells (`LAE`), fibroblasts (`FB`), and microvascular
endothelial cells (`MV`) — and then trained a neural network to predict, from the sequence
around a transcript's start and stop codons, whether NMD would destroy it.

This repository holds the code. The data it runs on is stored separately on Zenodo. Together
they rebuild the analyses the paper reports. What that rebuild reproduces exactly, and where it
differs, is recorded per claim rather than asserted here — see
[`docs/VERIFIED.md`](docs/VERIFIED.md).

---

## What you need

| what | why you need it |
|---|---|
| **This repository** | the analysis code |
| **The Zenodo archive** ([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) | the starting data — 12 data objects (7 files and 5 archives), 695 MB compressed, 4.6 GB once extracted |
| **The container image** *(in the same Zenodo archive)* | `nmd_1.2.sif`, 8.0 GB, carrying the exact R and Python versions the analyses ran under. Fetching it is optional: the `Dockerfile` here builds an equivalent environment (`docker build -t nmd:1.0 .`), restoring all 187 packages in `renv.lock` and both conda environments, and failing loudly if it cannot. See [`ENVIRONMENT.md`](ENVIRONMENT.md) |
| **`Isopair`** ([10.5281/zenodo.21536494](https://doi.org/10.5281/zenodo.21536494)) | an R package the analysis calls; `renv::restore()` fetches it from its GitHub record without a separate step |
| **`Isocall`** ([10.5281/zenodo.21536485](https://doi.org/10.5281/zenodo.21536485)) | the long-read pipeline that **produced** the deposited count matrices. Provenance, not a dependency — no shipped code calls it and it is not in `renv.lock` |
| **The model repository** `NMD_orf_model_v5_4ct` *(only for §5, and only if you are retraining)* | §5 is the one section that does not run from this repository alone — the training, evaluation and interpretation code lives there, and it needs a GPU. **You do not need it to reproduce §5's reported numbers**: the deposit's `model.zip` carries the trained weights, the per-isoform predictions and the full interpretability export, so every §5 figure and number rebuilds from the archive. Note that the public `NMD_orf_model_v5` is an **earlier, superseded model** and is not the one the paper reports |

The raw sequencing reads are in GEO under **GSE329233**. You do not need them. The Zenodo archive
starts one step later, from tables of how much of each transcript was seen in each sample, so you
can reproduce the analyses without reprocessing any sequencing data.

## Getting the data

Twelve data objects: seven you download directly, five that are zipped folders.

| file | what it contains |
|---|---|
| `nmd_lungcells_counts_4ct.csv` | how much of each transcript was seen in each sample, after quality filtering |
| `nmd_isocall_counts_4ct.csv` | the same counts over a slightly wider set of transcripts, before that filtering |
| `salmon_gene_counts_4ct.csv` | gene-level counts from short-read sequencing of the same samples, used as an independent check |
| `pheno_4ct.csv` | which sample is which — donor, cell type, and whether NMD was blocked |
| `nmd_isocall_4ct.gtf.gz` | the transcript models: where each transcript's exons sit on the genome |
| `encode_rbp_roster_vannostrand2020.csv` | the list of RNA-binding proteins used to test whether NMD targets are enriched for them |
| `gerstberger_2014_rbp_census.csv` | a second, independently compiled RNA-binding protein list, used as a check on the first |
| `sqanti.zip` | transcript classifications and predicted coding sequences (349 MB zipped, 4.0 GB open) |
| `isopair.zip` | the expressed-isoform universe and the four mashr differential-isoform tables the isoform-pair analysis starts from |
| `annotation.zip` | reference gene and transcript names from Ensembl release 115, and the gene-set table used for pathway enrichment |
| `model.zip` | the trained neural network and its interpretation outputs |
| `tan_2025_supplementary.zip` | published NMD target tables from another study, re-analysed here for comparison |

*(`4ct` in a filename just means the four cell types above.)*

Four more files sit alongside them in the archive: `MANIFEST.sha256`, which lets you check the
download; `UPLOAD_MANIFEST.sha256`, the checksums of the uploaded files as uploaded; the record's
own `README.md`, which describes the data in more detail than this page does; and the container
image `nmd_1.2.sif` described under **What you need** above.

Download them, unzip the five folders, and point the code at the result:

```bash
# from the repository root
ln -s /path/to/nmd_deposit_2026 data_deposit

# check the download — 58 lines, every one should say OK
cd data_deposit/source_data && shasum -a 256 -c ../MANIFEST.sha256 && cd -
```

The paths inside `MANIFEST.sha256` are relative to `source_data/`, which is why the check is run
from inside that folder rather than from the repository root.

That link is deliberately **not** stored in git. A stored one would carry the original author's
folder layout and then fail on your computer while looking as though it had worked.

## Running the analyses

[`REPRODUCTION.md`](REPRODUCTION.md) lists the steps in order with the command for each. Two
things are worth knowing before you start:

- **Some steps take hours, and training the neural network needs a graphics card.** The trained
  network is included in `model.zip` for exactly that reason — you do not have to retrain it, and
  everything after it runs on an ordinary laptop.
- **No script contains a hard-coded folder path.** They all read `config/paths.yml`. If a script
  cannot find its input, that file is the one to check.

Where a number we rebuilt differs from the one printed in the paper, the difference is written
down rather than quietly reconciled — [`docs/VERIFIED.md`](docs/VERIFIED.md) marks every such
claim `differs`.

## How the code is organised

This is the **minimal set**: the code needed to rebuild the paper's results, and nothing else. Our
own build and checking scripts, the working notes we kept while writing, and the intermediate
records are not here, because a reader does not need them to reproduce a number.

| folder | what is in it |
|---|---|
| `analysis/` | the main pipeline, in the order `REPRODUCTION.md` runs it |
| `figures/` | one folder per figure, holding the script for each panel and the script that assembles them |
| `R/`, `python/` | shared helper code — finding input files, consistent plot styling, and recording each reported number as it is computed |
| `config/` | where the input data lives (`paths.yml`) and which files are part of the analysis (`corpus.yml`) |
| `docs/` | one file only: [`docs/VERIFIED.md`](docs/VERIFIED.md), the per-claim record of what the rebuild reproduces. [`REPRODUCTION.md`](REPRODUCTION.md) sends you there |

## Citing this work

See [`CITATION.cff`](CITATION.cff). Cite the paper for the findings, the Zenodo source-data record
([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) for the data, and this package
for the code:

> **[10.5281/zenodo.21897099](https://doi.org/10.5281/zenodo.21897099)**

That is the concept DOI and always resolves to the latest version — cite it rather than a version
DOI, which keeps pointing at the files of the release that minted it.

## Licence

Released under the **MIT License** — see [`LICENSE`](LICENSE). `Isopair` and `Isocall` are
released separately with their own DOIs, and the source data and container image are archived on
Zenodo under that record's own terms (CC-BY-4.0).
