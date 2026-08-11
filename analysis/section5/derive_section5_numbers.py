#!/usr/bin/env python3
"""derive_section5_numbers.py -- every section 5 quantity, computed from the DEPOSITED model.

WHY THIS EXISTS. Section 5 is being replaced: the paper describes a 500nt, single (non-ensemble)
model, and the deposit now carries a 1000nt single model, seed 42, swapped 2026-08-05 (W321). So
every number in section 5 has to be re-derived rather than re-checked, and the manuscript edited to
match. This script is the producer for that re-derivation, so each value has a producer and a
population instead of living as prose.

CORRECTION 2026-08-06: an earlier version of this docstring said docs/VERIFIED.md's 5.2.1 values
were STALE. THAT WAS WRONG, and the guardian caught it. That field already records the
deposit-native reading -- "atg1000_stop1000 seed 42: AUC 0.9257, AUPRC 0.8175, n_test 10,522" --
together with the exact metrics JSON, and it already notes that the sweep re-selected 1000/1000. It
records FOUR readings of this claim deliberately, beside each other rather than over each other,
because overwriting one with another is how two quantities with the same name over different sets
get conflated. What is actually wrong there is ORDERING: the field LEADS with the superseded
post-clip 500/500 reading (job 8913382, 0.9191 / 0.8183, n=10,520), so a reader meets the dead
vintage first. That is a presentation defect, not a stale value, and the two should not be confused
-- which is the same distinction this script exists to enforce.

WHAT IT READS. Only $DEPOSIT (default data_deposit/source_data/model), i.e. the files a reader gets
from Zenodo. Nothing from ../NMD_orf_model_v5_4ct/results_4ct, which holds the SUPERSEDED 500nt
June vintage and whose metrics_atg500_stop500.json (AUC 0.9306 / AUPRC 0.8330 / n_test 10,131) is
the published number -- reading it here is how the old values would silently survive the swap.

CLAIM 5.3.3 COMES FROM A SECOND PRODUCER, NOT FROM HERE. model:10_export_stop_codon_freq_sf37.py
computes the stop-codon frequencies and their Fisher tests; this script READS its two output tables
and emits them. Deliberately not recomputed: two implementations of one quantity is the failure
this project is named for. That producer needs selected_orfs.tsv, which is an Explorer file, not a
deposit file -- fetched 2026-08-06 and verified by sha256 against the remote copy. Set
$STOP_CODON_DIR if the outputs are elsewhere; absent them, 5.3.3 is skipped and says so.

AGGREGATION. Branch shares use mean_of_abs: mean |SHAP| per branch, divided by the summed
per-branch means. That is the aggregation 11_kernel_shap_branches.py:562-566 implements, read out
of the producer rather than guessed -- and note that block computes over NMD ONLY, which is why the
published 60.7 / 28.8 / 10.5 matches the NMD-only universe and not the all-isoform one. Both are
emitted, per D78, so a value can never be read without the set it was computed over.

USAGE
    python analysis/section5/derive_section5_numbers.py            # emit
    python analysis/section5/derive_section5_numbers.py --print    # human table, no emission
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = next(p for p in HERE.parents if (p / "config" / "paths.yml").exists())
sys.path.insert(0, str(ROOT / "tools"))

TAG = "atg1000_stop1000_seed42"
DEPOSIT = Path(os.environ.get("DEPOSIT", ROOT / "data_deposit" / "source_data" / "model"))

# The universe strings are written once, here, because they are the part a reader needs and the
# part most easily lost: "58%" and "54.6%" are the same computation over different populations.
POP_TEST = ("held-out test_clean isoforms (chr 1/3/5/7, paralog-free), deposit-native 1000nt "
            "universe; training and validation excluded")
POP_ALL = ("all isoforms scored by the model, train and held-out pooled (D77/D74); "
           "interpretation only, never a performance number")
POP_NMD = "NMD susceptible isoforms only, train and held-out pooled -- the universe the PUBLISHED branch percentages were computed over"
POP_CTRL = "control (non-NMD-susceptible) isoforms only, train and held-out pooled"

# ARCHITECTURAL CONSTANTS ARE NOT COMPUTED OVER A POPULATION, AND SAYING SO IS THE FIX.
# The window size and the block count were emitted with the held-out test population and
# n = 10,522, which asserted they had been measured over the test set. They are properties of the
# model. Both seats then filed this as needing a fifth enumerated state and a contract change --
# and neither had read the eleven lines that settle it: `population` is FREE TEXT, validated for
# non-emptiness only (tools/claim_emit.py:191, R/claim_emit.R:48). There is no vocabulary to
# extend, so the honest string IS the fix. `n` is dropped for the same reason: there is no
# population to size.
#
# OPEN, and deliberately not resolved here: the field's own error text says "state what was
# counted", so a value that counted nothing is arguably not a claim_emit case at all. That is a
# design question about the contract rather than about these two values, it is not urgent, and the
# numbers are right and the field is now true. Raised by the guardian, 2026-08-06.
POP_ARCH = ("not applicable -- architectural constant, a property of the model rather than a "
            "quantity computed over any set of isoforms")


def _read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _require(path: Path) -> Path:
    """A missing input must raise. Constraint C3 of the redesign spec: the failure mode being
    guarded is a chunk that skips at exit 0 and reports nothing, which is how five section 5
    claims came to have no visible producer at all."""
    if not path.exists():
        raise SystemExit(
            f"MISSING DEPOSIT INPUT: {path}\n"
            f"  This script reads only the deposited model ({DEPOSIT}).\n"
            f"  It does NOT fall back to ../NMD_orf_model_v5_4ct/results_4ct -- that tree holds the\n"
            f"  superseded 500nt vintage, and falling back to it would reproduce the published\n"
            f"  numbers while appearing to have recomputed them."
        )
    return path


# ----------------------------------------------------------------------------- performance
def performance() -> dict:
    j = json.load(open(_require(DEPOSIT / f"metrics_{TAG}_test_clean.json")))
    preds = _read_tsv(_require(DEPOSIT / f"predictions_{TAG}_test_clean.tsv"))
    y = [int(float(r["label"])) for r in preds]
    s = [float(r["prob"]) for r in preds]
    n, npos = len(y), sum(y)
    if npos == 0 or npos == n:
        raise SystemExit("test split has one class only -- AUC undefined")

    order = sorted(range(n), key=lambda i: s[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j2 = i
        while j2 + 1 < n and s[order[j2 + 1]] == s[order[i]]:
            j2 += 1
        r = (i + j2) / 2 + 1
        for k in range(i, j2 + 1):
            ranks[order[k]] = r
        i = j2 + 1
    auc = (sum(ranks[i] for i in range(n) if y[i] == 1) - npos * (npos + 1) / 2) / (npos * (n - npos))

    tp = fp = 0
    prev_rec = 0.0
    ap = 0.0
    for i in sorted(range(n), key=lambda i: -s[i]):
        tp += y[i]
        fp += 1 - y[i]
        rec = tp / npos
        ap += (rec - prev_rec) * (tp / (tp + fp))
        prev_rec = rec

    # The recomputation is the check. If the shipped json and the predictions disagree, the pair is
    # inconsistent and neither can be quoted -- that is a defect, not a rounding difference.
    for name, mine, theirs in (("AUC", auc, j["auc"]), ("AUPRC", ap, j["auprc"])):
        if abs(mine - theirs) > 1e-3:
            raise SystemExit(f"{name} recomputed {mine:.4f} but metrics json says {theirs:.4f} -- "
                             f"the deposited json and the deposited predictions disagree")
    return dict(auc=auc, auprc=ap, n=n, n_nmd=npos, json=j)


# ----------------------------------------------------------------------------- branch shares
def branch_shares() -> dict:
    rows = _read_tsv(_require(DEPOSIT / f"kernel_shap_branch_{TAG}_all.tsv"))
    out = {}
    for name, sub in (("all", rows),
                      ("nmd", [r for r in rows if float(r["label"]) > 0.5]),
                      ("ctrl", [r for r in rows if float(r["label"]) < 0.5])):
        if not sub:
            raise SystemExit(f"branch universe '{name}' is empty -- an empty factor level renders "
                             f"without error and reports a share of nothing")
        m = {b: sum(abs(float(r[f"shap_{b}"])) for r in sub) / len(sub)
             for b in ("atg", "stop", "structural")}
        tot = sum(m.values())
        out[name] = dict(n=len(sub), mean=m,
                         pct={b: 100 * v / tot for b, v in m.items()},
                         stop_over_atg=m["stop"] / m["atg"])
    return out


# ----------------------------------------------------------------------------- structural features
def structural_features() -> dict:
    files = sorted(glob.glob(str(DEPOSIT / f"deepshap_summary_{TAG}_structural_run*.tsv")))
    if len(files) < 2:
        raise SystemExit(f"expected several DeepSHAP replicates, found {len(files)} -- the spread "
                         f"across replicates is the only handle on estimate noise")
    acc: dict[str, dict[str, list]] = {}
    for f in files:
        for r in _read_tsv(Path(f)):
            d = acc.setdefault(r["channel"], {"all": [], "nmd": [], "ctrl": []})
            d["all"].append(float(r["mean_abs_shap"]))
            d["nmd"].append(float(r["mean_abs_shap_nmd"]))
            d["ctrl"].append(float(r["mean_abs_shap_ctrl"]))
    ranked = sorted(acc, key=lambda k: -st.mean(acc[k]["all"]))
    return dict(n_replicates=len(files), acc=acc, ranked=ranked,
                lead_ratio=st.mean(acc[ranked[0]]["all"]) / st.mean(acc[ranked[1]]["all"]))


# ----------------------------------------------------------------------------- uORF attention
def uorf_attention() -> dict:
    rows = _read_tsv(_require(DEPOSIT / "uorf_attention_metrics.tsv"))
    g: dict[str, list] = {}
    for r in rows:
        g.setdefault(r["subgroup"], []).append(float(r["uorf_attention_frac"]))
    for k in ("NMD_PTC_neg_lost", "NMD_PTC_neg_retained", "NMD_PTC_pos"):
        if k not in g:
            raise SystemExit(f"subgroup {k} absent from uorf_attention_metrics.tsv")
    return {k: dict(n=len(v), median=st.median(v), mean=st.mean(v),
                    frac_nonzero=sum(1 for x in v if x > 0) / len(v))
            for k, v in g.items()}


# ----------------------------------------------------------------------------- main
# ----------------------------------------------------------------------------- stop codon usage
# Produced by model:10_export_stop_codon_freq_sf37.py, which needs selected_orfs.tsv -- an Explorer
# file, fetched 2026-08-06 and verified by sha256 against the remote copy (198,817 rows,
# 8f72aa93844c084f9ff32a5974e4b62f4f34271a5ec0f36b205ab5f28119e945). Read from its output rather
# than recomputed here: two implementations of one quantity is this project's standing failure.
# DEFAULTED TO THE DEPOSIT 2026-08-10. It used to default into the sibling model repo, so the
# clean room printed "5.3.3 NOT AVAILABLE" -- every section 5 number reproduced there except
# this one, for want of two files totalling 615 bytes. Both are now in the deposit. The env
# override stays for a run pointed at a fresh export.
STOP_DIR = Path(os.environ.get("STOP_CODON_DIR", DEPOSIT))


def stop_codon() -> dict | None:
    freq = STOP_DIR / "stop_codon_freq_by_class_sf37_all.tsv"
    test = STOP_DIR / "stop_codon_test_sf37_all.tsv"
    if not (freq.exists() and test.exists()):
        return None
    f = {(r["population"], r["stop_codon"]): r for r in _read_tsv(freq)}
    t = {r["stop_codon"]: r for r in _read_tsv(test)}
    return dict(freq=f, test=t)


# THE "pub" VALUES ARE THE CURRENT MANUSCRIPT, corrected 2026-08-10, and they were not.
#
# This script's comparison column is the thing a person reads to judge whether section 5
# reproduces. Four of its reference values were the SUPERSEDED 500nt model: AUPRC 0.83 (paper says
# 0.82), window 500nt (paper says 1000nt), lead ratio 2.153x (paper says 0.940 mean |SHAP|, ~3x),
# and "UGA 56% vs 48%, p<0.001" (paper says 58% vs 50%, OR 1.40, p 1.7e-42).
#
# The effect was to make a GOOD reconciliation look bad. Measured against the manuscript of
# 2026-07-17, the clean-room run (job 9057341) matches on every one of them: 0.9257/0.8175 -> 0.93
# and 0.82; 1000nt; 54.6/31.4/14.0 over 41,776; 0.9402 and 2.93x; and 5.3.3 exactly -- 58.1% vs
# 49.7%, OR 1.401, p 1.65e-42 against the paper's 58% vs 50%, OR 1.40, p 1.7e-42.
#
# The 500nt row is kept and RELABELLED rather than deleted: a superseded value beside its
# replacement is what stops it being re-derived, and this one has already leaked back twice.


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--print", dest="show", action="store_true",
                    help="human-readable table, emit nothing")
    args = ap.parse_args()

    perf = performance()
    br = branch_shares()
    sf = structural_features()
    ua = uorf_attention()
    sc = stop_codon()

    if args.show:
        j = perf["json"]
        print(f"MODEL {TAG}  windows {j['window_size_atg']}/{j['window_size_stop']}nt  "
              f"best_epoch {j['best_epoch']}  seed {j.get('member_seed')}")
        # pub 0.82, not 0.83. The manuscript reads "AUC of 0.93 and AUPRC of 0.82"; the emit call
        # below has always said published=0.82 and this print said 0.83, so the two disagreed inside
        # one script. It is the printed table a person reads to judge the reproduction, and against
        # the phantom 0.83 the deposited model's 0.8175 looked five times further off than it is.
        print(f"\n5.2.1  AUC {perf['auc']:.4f} (pub 0.93)   AUPRC {perf['auprc']:.4f} (pub 0.82)   "
              f"n_test {perf['n']:,} ({perf['n_nmd']:,} NMD)")
        print(f"       window {j['window_size_atg']}nt (pub 1000nt)")
        print("\n5.2.2 / 5.2.3  branch decomposition (mean_of_abs)")
        print(f"       {'universe':10s} {'n':>7s} {'ATG':>7s} {'stop':>7s} {'struct':>8s} {'stop/ATG':>9s}")
        for k, lbl in (("all", "all"), ("nmd", "NMD-only"), ("ctrl", "Control")):
            d = br[k]
            print(f"       {lbl:10s} {d['n']:7,d} {d['pct']['atg']:6.1f}% {d['pct']['stop']:6.1f}% "
                  f"{d['pct']['structural']:7.1f}% {d['stop_over_atg']:8.2f}x")
        print("       paper (Fig 5C, all 41,776):          14.0%   31.4%    54.6%")
        print("       SUPERSEDED 500nt model, NMD-only:    10.5%   28.8%    60.7%     2.74x")
        print(f"\n5.2.4  structural features, mean over {sf['n_replicates']} replicates")
        for k in sf["ranked"]:
            v = sf["acc"][k]
            print(f"       {k:20s} all {st.mean(v['all']):.4f}  NMD {st.mean(v['nmd']):.4f}")
        print(f"       lead ratio {sf['ranked'][0]} / {sf['ranked'][1]} = {sf['lead_ratio']:.2f}x "
              f"(paper 0.940 mean |SHAP|, ~3x the next feature)")
        print("\n5.4.3  uORF attention by subgroup")
        for k in ("NMD_PTC_neg_lost", "NMD_PTC_neg_retained", "NMD_PTC_pos", "non_NMD"):
            if k in ua:
                d = ua[k]
                print(f"       {k:22s} n={d['n']:6,d}  median {d['median']:.3f}  mean {d['mean']:.3f}")
        if sc is None:
            print("\n5.3.3  NOT AVAILABLE -- run model:10_export_stop_codon_freq_sf37.py "
                  "(needs selected_orfs.tsv from Explorer)")
        else:
            print("\n5.3.3  stop codon usage, priority (rank-0) ORF")
            for codon in ("UGA", "UAA", "UAG"):
                n = sc["freq"][("NMD", codon)]
                c = sc["freq"][("Control", codon)]
                t = sc["test"][codon]
                print(f"       {codon}  NMD {n['pct']:>5s}% ({n['n']:>6s}/{n['population_total']})   "
                      f"Control {c['pct']:>5s}% ({c['n']:>6s}/{c['population_total']})   "
                      f"OR {float(t['odds_ratio']):.3f}  p {float(t['p_value']):.2e}")
            print("       paper: UGA 58% vs 50%, OR 1.40, p 1.7e-42")
        return

    # claim_emit LEFT THIS REPO in ae89cc9 (D1020) and this call site was missed, because it is a
    # deferred import inside main() rather than a top-of-file one -- that commit's sweep concluded
    # "only TWO genuinely import them" and this is a third. The clean-room run of 2026-08-09 (job
    # 9045389) is where it surfaced: bare `python3 derive_section5_numbers.py` defaults to emit, so
    # the room got ModuleNotFoundError instead of a single section 5 number.
    #
    # The room does not get claim_tools and should not: emitting into the claim ledger is apparatus,
    # the ledger is not in the citable set, and D1020 moved these tools out deliberately. So the
    # chain runs --print, which computes and reports every quantity this path emits, and emit stays
    # the default here for the laptop-side ledger workflow. Say which mode to use rather than
    # dying on an import trace.
    try:
        from claim_emit import emit
    except ModuleNotFoundError as e:
        raise SystemExit(
            f"cannot emit: {e}. claim_emit lives in the claim_tools repository (D1020), not here.\n"
            "  To COMPUTE the section 5 numbers without it (this is what the rerun chain does):\n"
            "    python3 analysis/section5/derive_section5_numbers.py --print\n"
            "  To EMIT into the claim ledger, install claim_tools first."
        ) from e

    j = perf["json"]
    emit("5.2.1", "start/stop sequence window size (nt)", j["window_size_atg"], published=1000,
         population=POP_ARCH, restriction="none")
    emit("5.2.1", "test AUC", perf["auc"], published=0.93, n=perf["n"], population=POP_TEST,
         restriction="none")
    emit("5.2.1", "test AUPRC", perf["auprc"], published=0.82, n=perf["n"], population=POP_TEST,
         restriction="none")
    emit("5.2.1", "test set size", perf["n"], n=perf["n"], population=POP_TEST, restriction="none")
    emit("5.2.1", "NMD susceptible isoforms in test set", perf["n_nmd"], n=perf["n"],
         population=POP_TEST, restriction="none")

    emit("5.1.4", "number of input blocks", 3, published=3, population=POP_ARCH,
         restriction="none")

    # THE PUBLISHED COUNTERPART MOVED UNIVERSES ON 2026-08-08, and that is the whole point of
    # keeping it on the matching row. The paper USED to print 60.7 / 28.8 / 10.5 and 2.74, all
    # NMD-only and unlabelled. It now prints "55% of the predictive information" and "roughly
    # twice", both over all 41,776 scored isoforms (D77), and says so. So `published=` attaches to
    # the ALL rows and the NMD rows lose theirs -- they are still emitted, because D78 requires all
    # three universes recorded, but the paper no longer makes a claim about them.
    for key, pop in (("all", POP_ALL), ("nmd", POP_NMD), ("ctrl", POP_CTRL)):
        d = br[key]
        for branch in ("structural", "stop", "atg"):
            emit("5.2.2", f"{branch} branch share of total |SHAP| ({key})", d["pct"][branch],
                 published=55 if (key == "all" and branch == "structural") else None,
                 n=d["n"], population=pop, aggregation="mean_of_abs", restriction="none")
        emit("5.2.3", f"stop-window / start-window mean |SHAP| ratio ({key})", d["stop_over_atg"],
             published=2 if key == "all" else None, n=d["n"], population=pop,
             aggregation="mean_of_abs", restriction="none")

    # PUBLISHED VALUES GO ON THE MATCHING QUANTITY AND THE MATCHING UNIVERSE. This block used to
    # hang published=2.153 on the LEAD RATIO, which compared a computed ratio (2.93) against a
    # published MEAN -- two different quantities in one row, over two different universes. The
    # paper's legend carries both numbers and both are NMD-only: EJC mean |SHAP| 2.153, and
    # "~15x the next feature". So 2.153 belongs on an NMD-only mean and 15 on an NMD-only ratio,
    # and the all-isoform figures (D77) carry no published counterpart because the paper never
    # printed one. Caught by the guardian's normalization pass, 2026-08-06.
    lead, second = sf["ranked"][0], sf["ranked"][1]
    for ch in sf["ranked"]:
        emit("5.2.4", f"structural feature mean |SHAP| ({ch})", st.mean(sf["acc"][ch]["all"]),
             n=br["all"]["n"], population=POP_ALL, aggregation="mean_of_abs",
             sd_within=st.stdev(sf["acc"][ch]["all"]), replicate=sf["n_replicates"],
             restriction="none")
    emit("5.2.4", f"lead ratio, {lead} over {second}", sf["lead_ratio"],
         n=br["all"]["n"], population=POP_ALL, aggregation="mean_of_abs", restriction="none")

    # The NMD-only pair, which is what the published legend actually reports.
    nmd_rank = sorted(sf["acc"], key=lambda k: -st.mean(sf["acc"][k]["nmd"]))
    emit("5.2.4", f"structural feature mean |SHAP|, NMD only ({nmd_rank[0]})",
         st.mean(sf["acc"][nmd_rank[0]]["nmd"]), n=br["nmd"]["n"],
         population=POP_NMD, aggregation="mean_of_abs",
         sd_within=st.stdev(sf["acc"][nmd_rank[0]]["nmd"]), replicate=sf["n_replicates"],
         restriction="none")
    emit("5.2.4", f"lead ratio NMD only, {nmd_rank[0]} over {nmd_rank[1]}",
         st.mean(sf["acc"][nmd_rank[0]]["nmd"]) / st.mean(sf["acc"][nmd_rank[1]]["nmd"]),
         n=br["nmd"]["n"], population=POP_NMD, aggregation="mean_of_abs",
         restriction="none")

    for k, d in ua.items():
        emit("5.4.3", f"median uORF attention fraction ({k})", d["median"], n=d["n"],
             population=f"{k} isoforms, train and held-out pooled", restriction="none")

    if sc is None:
        print("NOT emitted: 5.3.3 stop-codon usage -- run model:10_export_stop_codon_freq_sf37.py",
              file=sys.stderr)
    else:
        # The published sentence quotes UGA only; all three codons are emitted so that quoting UGA
        # is visibly a choice rather than the only thing measured.
        for codon in ("UGA", "UAA", "UAG"):
            pub_n, pub_c = (58.0, 50.0) if codon == "UGA" else (None, None)
            for pop, pub in (("NMD", pub_n), ("Control", pub_c)):
                r = sc["freq"][(pop, codon)]
                emit("5.3.3", f"{codon} share of priority-ORF stop codons ({pop})",
                     float(r["pct"]), published=pub, n=int(r["population_total"]),
                     population=(f"{pop} isoforms with a valid rank-0 ORF stop codon, "
                                 f"train and held-out pooled"),
                     restriction="none")
            t = sc["test"][codon]
            emit("5.3.3", f"{codon} NMD-vs-Control odds ratio", float(t["odds_ratio"]),
                 published=1.40 if codon == "UGA" else None,
                 n=int(t["nmd_total"]) + int(t["control_total"]),
                 population="NMD vs Control, priority-ORF stop codon, Fisher exact",
                 restriction="none")
            emit("5.3.3", f"{codon} NMD-vs-Control Fisher p", float(t["p_value"]),
                 published=1.7e-42 if codon == "UGA" else None,
                 n=int(t["nmd_total"]) + int(t["control_total"]),
                 population="NMD vs Control, priority-ORF stop codon, Fisher exact",
                 restriction="none")

    print(f"emitted section 5 values from {DEPOSIT}", file=sys.stderr)


if __name__ == "__main__":
    main()
