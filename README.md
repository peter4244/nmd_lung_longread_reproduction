# Nonsense-mediated decay in human lung cells — analysis code

Code for the analyses and figures in *Long Read RNA Sequencing in Primary Lung Cell Types Reveals
Principles of Nonsense-Mediated Decay* (Leshem et al., 2026).

**The study in one paragraph.** Nonsense-mediated decay (NMD) is the cell's system for destroying
faulty messenger RNAs. We inhibited it in cultured human lung cells with an SMG1 inhibitor
(`Smg1i` in the data files, against a `DMSO` vehicle control), sequenced the RNA with long reads
so that whole transcripts could be read end to end, and compared what appeared when NMD was
inhibited against what appeared in untreated cells. That tells us which RNA molecules the cell had
been quietly destroying. We did this in four lung cell types — alveolar type 2 cells (`AT2`),
airway epithelial cells (`LAE`), fibroblasts (`FB`), and microvascular endothelial cells (`MV`) —
and then trained a neural network to predict, from the sequence around a transcript's start and
stop codons, whether NMD would destroy it.

This repository holds the code; the data it runs on is on Zenodo.

**Read in this order.** This page → [`ENVIRONMENT.md`](ENVIRONMENT.md) for the container →
[`REPRODUCTION.md`](REPRODUCTION.md) for the ordered run. Two more documents are worth knowing
about: the Zenodo record ships its **own `README.md`**, which describes the data and owns the
on-disk layout — it downloads with the record and is also shown on the record page — and
[`RESULTS_INDEX.md`](RESULTS_INDEX.md) maps reported numbers to the script that produces them,
their inputs, and the command to re-run them. It covers **59 of the manuscript's 149 sentences**
and states its own coverage — the rest carry no reported number or were excluded during triage.

---

## What you need

| What | Why |
|---|---|
| **This repository** | The analysis code. |
| **The Zenodo archive** ([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) | The starting data — 12 objects, 695 MB compressed, 4.65 GB extracted. |
| **The container** *(same archive)* | `nmd_1.3.sif` (8.1 GB), or the same environment as a Docker archive (8.8 GB). Optional — the `Dockerfile` here builds an equivalent image. See [`ENVIRONMENT.md`](ENVIRONMENT.md). |
| **`Isopair`** ([10.5281/zenodo.21536494](https://doi.org/10.5281/zenodo.21536494)) | An R package the analysis calls. `renv::restore()` fetches it; no separate step. |
| **`Isocall_v1`** ([10.5281/zenodo.21536485](https://doi.org/10.5281/zenodo.21536485)) | The pipeline that produced the deposited count matrices. Provenance only — no shipped code calls it, and you do not need it to reproduce anything. |
| **`NMD_orf_model_v5_4ct`** ([10.5281/zenodo.21536501](https://doi.org/10.5281/zenodo.21536501)) | **Required, and not only for §5.** `export_rds.R` lives there and is the last step of §5.3, which runs *before* the §4 report — and it needs a batch allocation of roughly 96 GB, not a laptop. Retraining §5 additionally needs a GPU; §5's *reported numbers* rebuild from `model.zip` without one. |

> **Use `NMD_orf_model_v5_4ct`, not `NMD_orf_model_v5`.** The latter is an earlier, superseded
> model and is not the one the paper reports.

## Getting the data

Download the record ([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) and
follow **its own `README.md`**, which owns the on-disk layout: what each file contains, the
arrange-and-verify commands, and the checksum manifests. The record is open access — no token or
account needed.

One hazard is worth repeating here. The code resolves its deposit root to `<deposit>/source_data`,
so **extracting flat leaves every path pointing into a directory that does not exist.** The wiring
check in [`REPRODUCTION.md`](REPRODUCTION.md) §3 catches this and warns `DOES NOT RESOLVE`, but it
warns rather than stops — so read its output instead of assuming a clean run.

[`REPRODUCTION.md`](REPRODUCTION.md) picks up from there, starting with the symlink that points
this repository at your copy of the data.

## Running the analyses

[`REPRODUCTION.md`](REPRODUCTION.md) gives the steps in order with the command for each. Two
things to know first:

- **Some steps take hours, and two need a cluster.** Training the network needs a GPU — the trained
  network ships in `model.zip` so you do not have to. But `export_rds.R` in §5.3 needs a batch
  allocation of roughly 96 GB and is killed on a login node, so the chain is not laptop-only.
- **No script hard-codes a folder path.** They all read `config/paths.yml`. If a script cannot
  find its input, check that file.

## How the code is organized

This is the minimal set: the code needed to rebuild the paper's results, and nothing else.

| Folder | Contents |
|---|---|
| `analysis/` | The main pipeline, in the order `REPRODUCTION.md` runs it |
| `figures/` | One folder per figure — the script for each panel and the script that assembles them |
| `R/`, `python/` | Shared helpers: input-file resolution, plot styling, and recording each reported number as it is computed |
| `config/` | Where the data lives (`paths.yml`) |
| `RESULTS_INDEX.md` | Reported numbers mapped to their producing script, required inputs, and re-run command — 59 of the manuscript's 149 sentences, coverage stated in the file |

## Citing

See [`CITATION.cff`](CITATION.cff). Cite the paper for the findings, the Zenodo source-data record
([10.5281/zenodo.21544336](https://doi.org/10.5281/zenodo.21544336)) for the data, and this
package for the code:

> **[10.5281/zenodo.21897099](https://doi.org/10.5281/zenodo.21897099)**

That is the concept DOI and always resolves to the latest version — cite it rather than a version
DOI, which keeps pointing at the files of the release that minted it.

## License

MIT — see [`LICENSE`](LICENSE). `Isopair` and `Isocall_v1` are released separately under their own
DOIs, and the archived data and container images are on Zenodo under that record's terms
(CC-BY-4.0).
