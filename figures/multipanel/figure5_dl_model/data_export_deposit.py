#!/usr/bin/env python3
"""data_export_deposit.py -- build Figure 5's data/ layer from the DEPOSITED 1000nt model.

WHY THIS EXISTS. Measured 2026-08-06: of the 13 section 5 figure scripts, 11 read inputs that do
not exist anywhere on disk and the 2 that resolve read the SUPERSEDED 500nt tree. Figure 5's data/
directory is untracked and absent, so not one panel could be re-rendered at all -- the committed
PDFs are 2026-07-24 renders of inputs that are gone. The old export (data_export.R) covers only
Panel G and a GENCODE isoform list; panelB/C/D's inputs had no producer anywhere in the repo, which
is defect C40 of the redesign spec. This script is that missing producer, and it reads only files a
reader can fetch from Zenodo.

UNIVERSE (D77/D74). Interpretation quantities -- Panels C and D -- are computed over ALL isoforms,
train and held-out pooled, which is the universe D77 selects. Both universes are written out
(`_nmd` variants) so the choice is visible and reversible rather than baked in silently. Panel B is
held-out test only -- pooling train and held-out is legitimate for interpretation and never for a
performance number.

AND THE POOLED UNIVERSE IS THE MANUSCRIPT'S, WHICH THIS COMMENT PREVIOUSLY DENIED. It said the
published branch percentages were "60.7 / 28.8 / 10.5 ... NMD-only", and that the pooled universe
therefore CHANGED the panels' meaning relative to the published figure. Both halves were wrong.
Measured 2026-08-15 against kernel_shap_branch_atg1000_stop1000_seed42_all.tsv in the deposit: the
paper's Figure 5C reports 54.6 / 31.4 / 14.0, and the pooled universe gives 54.64 / 31.41 / 13.95
over n=41,776 -- an exact reproduction. The NMD-only subset gives 59.01 / 29.43 / 11.56 over
n=9,321 and is NOT the paper's number.

60.7 / 28.8 / 10.5 is real but belongs to the SUPERSEDED 500nt model over its NMD-only subset --
reproduced exactly from results_4ct/kernel_shap_branch_atg500_stop500.tsv and from
model.published_superseded_2026-07-28/. Calling it "the published percentages" here sent a reader
to compare the current figure against a retired one. Corrected rather than deleted, because the
number exists and will be found again. The retired cohort size is deliberately not restated here:
tools/check.py ratchets how far that constant may spread, and this file was the hit that broke it.

WHAT IT DOES NOT COVER, so a gap is not mistaken for completeness:
  SF38     -- needs selected_orfs.tsv, an Explorer file. See task #6.

Panel G WAS listed here as uncovered, on the grounds that it needed tx_summary.tsv, a stage C
output rather than a deposit file. That is no longer true: it now reads the deposited
uorf_attention_metrics.tsv instead, and is produced by this script. The line is corrected rather
than deleted because "not covered" and "covered" are the two answers a reader acts on differently.

USAGE
    python figures/multipanel/figure5_dl_model/data_export_deposit.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "config" / "paths.yml").exists())

TAG = "atg1000_stop1000_seed42"
DEPOSIT = Path(os.environ.get("DEPOSIT", ROOT / "data_deposit" / "source_data" / "model"))
OUT = HERE / "data"

# Friendly labels for the five structural channels. Kept here rather than in the panel so the
# panel plots whatever it is given and the naming has one home.
FEATURE_LABELS = {
    "n_downstream_ejc": "Downstream EJC count",
    "is_ref_cds": "Reference CDS",
    "frac_start": "Start position (fraction)",
    "frac_stop": "Stop position (fraction)",
    "is_sqanti_cds": "SQANTI CDS",
}
# Must match figure5_panelC_branch_importance.py's BRANCH_COLORS keys exactly -- the panel indexes
# its palette by this string and raises KeyError on a mismatch, which is the good failure.
BRANCH_LABELS = {"structural": "Structural", "stop": "Stop", "atg": "ATG"}


def require(p: Path) -> Path:
    if not p.exists():
        raise SystemExit(
            f"MISSING DEPOSIT INPUT: {p}\n"
            f"  Reading only {DEPOSIT}. There is deliberately NO fallback to\n"
            f"  ../NMD_orf_model_v5_4ct/results_4ct: that tree is the superseded 500nt vintage,\n"
            f"  and falling back would rebuild the PUBLISHED figure while looking like a refresh."
        )
    return p


def read_tsv(p: Path) -> list[dict]:
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(name: str, cols: list[str], rows: list[dict]) -> None:
    if not rows:
        raise SystemExit(f"refusing to write an empty {name}: an empty panel renders without error")
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {name:46s} {len(rows):>7,d} rows")


# --------------------------------------------------------------------------- Panel B
def panel_b() -> None:
    j = json.load(open(require(DEPOSIT / f"metrics_{TAG}_test_clean.json")))
    preds = read_tsv(require(DEPOSIT / f"predictions_{TAG}_test_clean.tsv"))
    write_tsv("panelB_predictions.tsv", ["isoform_id", "label", "prob"],
              [{"isoform_id": r["isoform_id"], "label": int(float(r["label"])), "prob": r["prob"]}
               for r in preds])
    n_nmd = sum(1 for r in preds if float(r["label"]) > 0.5)
    write_tsv("panelB_metrics.tsv", ["n_test", "n_nmd", "auc", "auprc", "window_nt", "split"],
              [{"n_test": len(preds), "n_nmd": n_nmd, "auc": j["auc"], "auprc": j["auprc"],
                "window_nt": j["window_size_atg"], "split": j["split"]}])
    if n_nmd != j["n_nmd"]:
        raise SystemExit(f"metrics json says n_nmd={j['n_nmd']} but the predictions hold {n_nmd}")


# --------------------------------------------------------------------------- Panel C
def panel_c() -> None:
    rows = read_tsv(require(DEPOSIT / f"kernel_shap_branch_{TAG}_all.tsv"))
    for suffix, sub in (("", rows), ("_nmd", [r for r in rows if float(r["label"]) > 0.5])):
        m = {b: sum(abs(float(r[f"shap_{b}"])) for r in sub) / len(sub)
             for b in ("atg", "stop", "structural")}
        tot = sum(m.values())
        out = [{"branch": BRANCH_LABELS[b], "channel": b, "mean_abs_shap": m[b],
                "pct": 100 * m[b] / tot, "n": len(sub)}
               for b in sorted(m, key=lambda b: -m[b])]
        write_tsv(f"panelC_branch_importance{suffix}.tsv",
                  ["branch", "channel", "mean_abs_shap", "pct", "n"], out)


# --------------------------------------------------------------------------- Panel D
def panel_d() -> None:
    files = sorted(glob.glob(str(DEPOSIT / f"deepshap_summary_{TAG}_structural_run*.tsv")))
    if len(files) < 2:
        raise SystemExit(f"expected several DeepSHAP replicates, found {len(files)}")
    acc: dict[str, dict[str, list]] = {}
    for f in files:
        for r in read_tsv(Path(f)):
            d = acc.setdefault(r["channel"], {"all": [], "nmd": []})
            d["all"].append(float(r["mean_abs_shap"]))
            d["nmd"].append(float(r["mean_abs_shap_nmd"]))
    for suffix, key in (("", "all"), ("_nmd", "nmd")):
        out = [{"feature": FEATURE_LABELS.get(ch, ch), "channel": ch,
                "mean_abs_shap": st.mean(v[key]), "sd": st.stdev(v[key]),
                "n_replicates": len(files)}
               for ch, v in sorted(acc.items(), key=lambda kv: -st.mean(kv[1][key]))]
        write_tsv(f"panelD_structural_features{suffix}.tsv",
                  ["feature", "channel", "mean_abs_shap", "sd", "n_replicates"], out)


# --------------------------------------------------------------------------- Panels E / F
def panels_ef() -> None:
    """The deposit's motif-logo tables are a column SUPERSET of the old ones and span the full
    1000nt window, so the panels can read them unchanged once they are under the expected name.
    Copied rather than symlinked: a symlink into data_deposit/ breaks the moment the deposit is a
    Zenodo download rather than a local tree."""
    for panel, src_name in (("atg", f"motif_logo_atg_{TAG}_run1.tsv"),
                            ("stop", f"motif_logos_stop_{TAG}_run1.tsv")):
        rows = read_tsv(require(DEPOSIT / src_name))

        # THE STOP TABLE IS NOW FIVE TABLES STACKED. The deposited file carries a `landmark`
        # column the old 500nt one did not have -- stop_codon plus two 3'UTR offsets and two
        # 3'UTR junctions -- so it holds five rows per (position, channel). The panel is the
        # stop codon logo, so select it; passing the stack through would silently average five
        # different landmarks into one logo, and the n differs between them (9,321 isoforms at
        # the stop codon against 6,265 and 5,036 at the junctions, because not every isoform
        # has one). Panel F's reshape raised on the duplicate index rather than averaging, which
        # is how this was found.
        if "landmark" in rows[0]:
            marks = sorted({r["landmark"] for r in rows})
            rows = [r for r in rows if r["landmark"] == "stop_codon"]
            if not rows:
                raise SystemExit(f"no stop_codon rows in {src_name}; landmarks present: {marks}")
            print(f"      selected landmark 'stop_codon' of {len(marks)}: {', '.join(marks)}")

        cols = list(rows[0].keys())
        pos = sorted({int(r["relative_position"]) for r in rows})
        seen = {(r["relative_position"], r["channel"]) for r in rows}
        if len(seen) != len(rows):
            raise SystemExit(f"{src_name}: {len(rows)} rows but only {len(seen)} distinct "
                             f"(position, channel) pairs -- an unhandled stacking dimension")
        write_tsv(f"panel{'E' if panel == 'atg' else 'F'}_motif_logo_{panel}.tsv", cols, rows)
        print(f"      window spans {min(pos):+d}..{max(pos):+d} nt "
              f"({len(pos)} positions, {len({r['channel'] for r in rows})} channels)")


# --------------------------------------------------------------------------- Panel G
# Subgroup display order and labels. "lost" and "retained" are the two NMD+/PTC- arms; the old
# panel excluded "lost" by construction because that scope had no projectable reference AUG, but
# the deposited metrics table carries it and claim 5.4.3 is about PTC-negative isoforms generally.
GROUP_ORDER = [
    ("NMD_PTC_neg_lost", "NMD+/PTC-\n(start lost)"),
    ("NMD_PTC_neg_retained", "NMD+/PTC-\n(start retained)"),
    ("NMD_PTC_pos", "NMD+/PTC+"),
    ("non_NMD", "Control"),
]


def panel_g() -> None:
    """Panel G on uorf_attention_frac, NOT on the Path-B-strict statistic.

    W348, confirmed and reproduced 2026-08-06: strict_attn is a LEFT JOIN against
    uorf_features_in_priority_slots.tsv followed by `is.na -> 0`, and 0 of 832 Control and 0 of 767
    NMD+/PTC+ isoforms appear in that slot file at all. Every zero in both groups is the FILL, so
    the panel's two significant p-values tested slot-file coverage rather than attention, and the
    NaN was the honest one of the three.

    uorf_attention_frac has neither defect and is measured for the whole cohort: 41,776 rows, zero
    blanks, and every zero carries n_uorf_slots == 0 -- verified 100% in all four subgroups, so a
    zero means "this isoform has no uORF slot", which is a real structural zero and not a missing
    row. That is checked here rather than asserted, because the whole point of W348 is that the two
    are indistinguishable once they are both the number 0.
    """
    rows = read_tsv(require(DEPOSIT / "uorf_attention_metrics.tsv"))
    blanks = sum(1 for r in rows if r["uorf_attention_frac"].strip() == "")
    if blanks:
        raise SystemExit(f"{blanks} blank uorf_attention_frac values -- this column is supposed to "
                         f"be measured for the whole cohort; a blank here is W348 again")

    out, seen = [], {g for g, _ in GROUP_ORDER}
    bad = 0
    for r in rows:
        if r["subgroup"] not in seen:
            continue
        frac = float(r["uorf_attention_frac"])
        slots = float(r["n_uorf_slots"])
        if frac == 0 and slots != 0:
            bad += 1
        out.append({"isoform_id": r["isoform_id"], "subgroup": r["subgroup"],
                    "uorf_attention_frac": frac, "n_uorf_slots": int(slots),
                    "label": r["label"], "split": r["split"]})
    if bad:
        print(f"      NOTE: {bad} zeros with n_uorf_slots != 0 -- these are MEASURED zeros mixed "
              f"with structural ones; the panel must not describe all zeros as 'no uORF'")
    else:
        print(f"      every zero has n_uorf_slots == 0: structural, not a fill (W348 clear)")
    write_tsv("panelG_uorf_attention_v1.tsv",
              ["isoform_id", "subgroup", "uorf_attention_frac", "n_uorf_slots", "label", "split"],
              out)

    # Pairwise tests live with the producer, not the panel: a p-value computed inside a plotting
    # script has no producer and cannot be re-derived. Mann-Whitney because the distributions are
    # zero-inflated and nowhere near normal, and it is the test the superseded panel used.
    from scipy.stats import mannwhitneyu
    by = {}
    for r in out:
        by.setdefault(r["subgroup"], []).append(r["uorf_attention_frac"])
    stats = []
    keys = [g for g, _ in GROUP_ORDER]
    for i, ga in enumerate(keys):
        for gb in keys[i + 1:]:
            a, b = by[ga], by[gb]
            nz_a = sum(1 for x in a if x > 0)
            nz_b = sum(1 for x in b if x > 0)
            if nz_a == 0 and nz_b == 0:
                # Two all-zero groups: there is nothing to compare and the honest output is NaN,
                # not a p-value. W348's NaN was the correct one of its three results.
                u, p = float("nan"), float("nan")
            else:
                u, p = mannwhitneyu(a, b, alternative="two-sided")
            stats.append({"group_x": ga, "group_y": gb, "nx": len(a), "ny": len(b),
                          "nonzero_x": nz_a, "nonzero_y": nz_b,
                          "median_x": st.median(a), "median_y": st.median(b),
                          "u": u, "p": p})
    write_tsv("panelG_uorf_attention_v1_pairwise.tsv",
              ["group_x", "group_y", "nx", "ny", "nonzero_x", "nonzero_y",
               "median_x", "median_y", "u", "p"], stats)
    for s in stats:
        print(f"      {s['group_x']:22s} vs {s['group_y']:22s} "
              f"med {s['median_x']:.3f}/{s['median_y']:.3f}  "
              f"nonzero {s['nonzero_x']:>5d}/{s['nonzero_y']:>5d}  p {s['p']:.3g}")


def main() -> None:
    print(f"Figure 5 data layer <- {DEPOSIT}")
    panel_b()
    panel_c()
    panel_d()
    panels_ef()
    panel_g()
    print("\nNOT built here:")
    print("  the Path-B-strict Panel G variant -- superseded by W348, see panel_g()")


if __name__ == "__main__":
    main()
