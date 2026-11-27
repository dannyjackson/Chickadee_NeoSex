#!/usr/bin/env python3
"""
8b_dxy_summary.py: turn pixy windows from 8a into the Z-vs-autosome dXY test.

    python3 8b_dxy_summary.py ${PROJDIR}/analyses/dxy/pixy --out ${PROJDIR}/analyses/dxy/dxy_Z

Reads <pixydir>/Z_{pi,dxy}.txt and, if present, <pixydir>/autosomes_{pi,dxy}.txt.

What it does
  * Labels each Z window by compartment, using the PhaseWY / Genespace boundaries
    (BCCH reference coordinates). Windows crossing a boundary are labelled "boundary"
    and left out of the compartment totals.
  * Totals are sum(count_diffs) / sum(count_comparisons), never a mean of window averages.
  * Mutation-rate scaling: R = mean dXY between each ingroup species and the outgroup (BOCH),
    computed per compartment. dXY / R removes the faster Z rate (male mutation bias).
  * Net divergence Da = dXY - (pi1 + pi2) / 2, roughly 2*mu*T once ancestral diversity is removed.
    Years = Da / (2 * mu_c) * g, where mu_c = mu * R_c / R_autosomes and g is the pair's mean
    generation length.
  * 95% CIs from a block bootstrap (1 Mb blocks by default).
  * The test: scaled dXY between BCCH and MOCH on the Z, divided by scaled dXY between BCCH and
    CACH on the autosomes.
        model 1 (BCCH-MOCH introgression of the Z after the BCCH-CACH split): ratio < 1
        model 2 (BCCH-MOCH species tree, autosomes swamped by BCCH-CACH gene flow): ratio > 1
    The Z's smaller ancestral Ne (3/4) pulls the dXY ratio down slightly under either model;
    the Da ratio is less affected by this.

Outputs: <out>.windows.tsv, <out>.summary.tsv, <out>.test.tsv, <out>.png
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

# BCCH-reference coordinates of the Z (CM059970.1 = NC_081289.1), from PhaseWY / Genespace.
Z_REGIONS = [
    (0,           14_782_500,  "neo_PAR"),            # PAR in BCCH and MOCH; autosomal in CACH
    (14_782_500,  25_440_000,  "MOCH_stratum"),       # PAR in BCCH, non-PAR in MOCH
    (25_440_000,  62_170_000,  "neo_PAR"),
    (62_170_000,  63_070_000,  "MOCH_stratum"),
    (63_070_000,  65_197_050,  "neo_nonPAR"),         # non-PAR in both BCCH and MOCH
    (65_197_050,  65_501_052,  "uncertain"),
    (65_501_052,  70_940_000,  "anc_nonPAR_CACHgap"), # ancestral non-PAR; not placed on CACH Z
    (70_940_000,  146_352_500, "anc_nonPAR"),
    (146_352_500, 10**9,       "anc_PAR"),            # PAR in MOCH and CACH; none found in BCCH
]
LABEL = {
    "autosomes": "Autosomes", "autosomes_lowpi": "Autosomes (low π)",
    "neo_PAR": "Neo-PAR", "MOCH_stratum": "MOCH-only stratum", "neo_nonPAR": "Neo non-PAR",
    "uncertain": "Uncertain", "anc_nonPAR_CACHgap": "Anc. non-PAR (CACH gap)",
    "anc_nonPAR": "Ancestral non-PAR", "anc_PAR": "Ancestral PAR",
}
ORDER = ["BCCH", "CACH", "MOCH", "BOCH"]
PAIRCOL = {"BCCH_CACH": "#6a4ce2", "BCCH_MOCH": "#008a17", "CACH_MOCH": "#ec4647"}  # as in 6e
PAIRNAME = {"BCCH_CACH": "Black-capped vs Carolina", "BCCH_MOCH": "Black-capped vs Mountain",
            "CACH_MOCH": "Carolina vs Mountain"}
TEXT2 = "#666666"
NONPAR = {"neo_nonPAR", "anc_nonPAR", "anc_nonPAR_CACHgap"}   # compartments where the test applies
MIN_BLOCKS = 5


def parse_gen(s):
    return {k.upper(): float(v) for k, v in (x.split("=") for x in s.split(","))}


def pair_key(a, b):
    a, b = sorted((a, b), key=lambda x: ORDER.index(x) if x in ORDER else 99)
    return f"{a}_{b}"


def read_pixy(pixydir, prefix):
    """Wide table: one row per window, count columns for each pi pop and each dXY pair."""
    fd, fp = (os.path.join(pixydir, f"{prefix}_{s}.txt") for s in ("dxy", "pi"))
    if not (os.path.exists(fd) and os.path.exists(fp)):
        return None
    key = ["chromosome", "window_pos_1", "window_pos_2"]
    d = pd.read_csv(fd, sep="\t", na_values=["NA"])
    d["pair"] = [pair_key(a, b) for a, b in zip(d.pop1, d.pop2)]
    d = d.pivot_table(index=key, columns="pair", values=["count_diffs", "count_comparisons", "no_sites"],
                      aggfunc="sum")
    d.columns = [f"{s}:{p}" for s, p in d.columns]
    p = pd.read_csv(fp, sep="\t", na_values=["NA"])
    p = p.pivot_table(index=key, columns="pop", values=["count_diffs", "count_comparisons"], aggfunc="sum")
    p.columns = [f"pi_{s}:{q}" for s, q in p.columns]
    w = d.join(p, how="outer").fillna(0).reset_index()
    return w.rename(columns={"chromosome": "chrom", "window_pos_1": "start", "window_pos_2": "end"})


def z_compartment(s, e):
    for a, b, lab in Z_REGIONS:
        if s >= a and e <= b + 1:          # pixy windows are 1-based, inclusive
            return lab
    return "boundary"


def stats(S, pairs, pops, out, R_auto=None, mu=None, gen=None):
    """S: dict column -> array of summed counts (bootstrap reps on axis 0). Returns dict of arrays."""
    with np.errstate(divide="ignore", invalid="ignore"):
        dxy = {p: S[f"count_diffs:{p}"] / S[f"count_comparisons:{p}"] for p in pairs}
        pi = {q: S[f"pi_count_diffs:{q}"] / S[f"pi_count_comparisons:{q}"] for q in pops}
        ing = [q for q in pops if q != out]
        R = np.mean([dxy[pair_key(q, out)] for q in ing], axis=0)
        res = {"R": R}
        for p in pairs:
            if out in p.split("_"):
                continue
            a, b = p.split("_")
            da = dxy[p] - (pi[a] + pi[b]) / 2
            res[f"dxy:{p}"] = dxy[p]
            res[f"dxy_n:{p}"] = dxy[p] / R
            res[f"da_n:{p}"] = da / R
            if R_auto is not None:
                mu_c = mu * R / R_auto
                res[f"years:{p}"] = da / (2 * mu_c) * np.mean([gen[a], gen[b]])
        return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pixydir")
    ap.add_argument("--out", default="dxy_Z")
    ap.add_argument("--outgroup", default="BOCH")
    ap.add_argument("--mu", type=float, default=4.6e-9, help="autosomal mutation rate per site per generation")
    ap.add_argument("--gen", default="bcch=3.0,cach=2.3,moch=2.6")
    ap.add_argument("--min-sites", type=int, default=10000, help="minimum callable sites per window")
    ap.add_argument("--block", type=float, default=1e6, help="bootstrap block size (bp)")
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--lowpi-quantile", type=float, default=0.25,
                    help="autosomal windows with mean ingroup pi at or below this quantile = low-recombination proxy")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    gen = parse_gen(a.gen)
    out = a.outgroup.upper()
    rng = np.random.default_rng(a.seed)

    z = read_pixy(a.pixydir, "Z")
    if z is None:
        sys.exit(f"No Z_dxy.txt / Z_pi.txt in {a.pixydir}")
    z["compartment"] = [z_compartment(s, e) for s, e in zip(z.start, z.end)]
    auto = read_pixy(a.pixydir, "autosomes")

    pairs = sorted({c.split(":")[1] for c in z.columns if c.startswith("count_diffs:")})
    pops = sorted({c.split(":")[1] for c in z.columns if c.startswith("pi_count_diffs:")},
                  key=lambda x: ORDER.index(x) if x in ORDER else 99)
    if out not in pops:
        sys.exit(f"Outgroup {out} not among pixy populations {pops}")
    ingroup_pairs = [p for p in pairs if out not in p.split("_")]
    cols = [c for c in z.columns if c.startswith(("count_", "pi_count_"))]

    frames = [z]
    if auto is not None:
        ing = [q for q in pops if q != out]
        with np.errstate(divide="ignore", invalid="ignore"):
            mpi = np.mean([auto[f"pi_count_diffs:{q}"] / auto[f"pi_count_comparisons:{q}"] for q in ing], axis=0)
        auto["compartment"] = "autosomes"
        thr = np.nanquantile(mpi, a.lowpi_quantile)
        low = auto[mpi <= thr].copy()
        low["compartment"] = "autosomes_lowpi"
        frames += [auto, low]
    W = pd.concat(frames, ignore_index=True)

    # drop thin windows (any pair below the minimum number of callable sites)
    site_cols = [f"no_sites:{p}" for p in pairs]
    W["min_sites"] = W[site_cols].min(axis=1)
    W = W[W.min_sites >= a.min_sites].copy()
    W["block"] = W.chrom + ":" + (W.start // a.block).astype(int).astype(str)

    # per-window scaled values for plotting
    with np.errstate(divide="ignore", invalid="ignore"):
        Rw = np.mean([W[f"count_diffs:{pair_key(q, out)}"] / W[f"count_comparisons:{pair_key(q, out)}"]
                      for q in pops if q != out], axis=0)
        for p in ingroup_pairs:
            W[f"dxy:{p}"] = W[f"count_diffs:{p}"] / W[f"count_comparisons:{p}"]
            W[f"dxy_n:{p}"] = W[f"dxy:{p}"] / Rw
    W.to_csv(f"{a.out}.windows.tsv", sep="\t", index=False)

    comps = [c for c in ["autosomes", "autosomes_lowpi", "neo_PAR", "MOCH_stratum", "neo_nonPAR",
                         "anc_nonPAR_CACHgap", "anc_nonPAR", "anc_PAR"] if c in set(W.compartment)]

    # point estimates and block-bootstrap draws
    def summed(sub, boot):
        # resample 1 Mb blocks; a compartment too small for that (the 2 Mb neo non-PAR) falls back
        # to single windows, which ignores autocorrelation, so its CI is too narrow (flagged below)
        unit = "block" if sub.block.nunique() >= MIN_BLOCKS else "start"
        B = sub.groupby(unit)[cols].sum()
        M = B.to_numpy(float)
        pt = dict(zip(cols, M.sum(axis=0)))
        idx = rng.integers(0, len(M), size=(boot, len(M)))
        bs = M[idx].sum(axis=1)
        return pt, {c: bs[:, i] for i, c in enumerate(cols)}, (len(B) if unit == "block" else -len(B))

    R_auto_pt = R_auto_bs = None
    if "autosomes" in comps:
        pt, bs, _ = summed(W[W.compartment == "autosomes"], a.boot)
        R_auto_pt = stats(pt, pairs, pops, out)["R"]
        R_auto_bs = stats(bs, pairs, pops, out)["R"]

    rows, draws = [], {}
    for c in comps:
        sub = W[W.compartment == c]
        pt, bs, nb = summed(sub, a.boot)
        s_pt = stats(pt, pairs, pops, out, R_auto_pt, a.mu, gen)
        s_bs = stats(bs, pairs, pops, out, R_auto_bs, a.mu, gen)
        draws[c] = s_bs
        for p in ingroup_pairs:
            r = {"compartment": c, "pair": p, "n_windows": len(sub),
                 "ci_unit": "1Mb blocks" if nb > 0 else "windows (CI too narrow)", "n_units": abs(nb),
                 "R_outgroup_dxy": s_pt["R"]}
            for k in ("dxy", "dxy_n", "da_n", "years"):
                key = f"{k}:{p}"
                if key in s_pt:
                    lo, hi = np.nanpercentile(s_bs[key], [2.5, 97.5])
                    r[k], r[f"{k}_lo"], r[f"{k}_hi"] = s_pt[key], lo, hi
            rows.append(r)
        draws[c]["_pt"] = s_pt
        draws[c]["_nblocks"] = nb
    summ = pd.DataFrame(rows)
    summ.to_csv(f"{a.out}.summary.tsv", sep="\t", index=False, float_format="%.6g")

    # the test: Z BCCH-MOCH vs autosomal BCCH-CACH, both outgroup-scaled
    tests = []
    if "autosomes" in comps:
        for zc in [c for c in comps if not c.startswith("autosomes")]:
            for ac in [c for c in comps if c.startswith("autosomes")]:
                for k in ("dxy_n", "da_n"):
                    num, den = f"{k}:BCCH_MOCH", f"{k}:BCCH_CACH"
                    pt = draws[zc]["_pt"][num] / draws[ac]["_pt"][den]
                    bs = draws[zc][num] / draws[ac][den]
                    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
                    # the test only applies where BCCH and MOCH do not recombine their Z with the W;
                    # in recombining compartments (PARs) the ratio is shown for reference only
                    if zc not in NONPAR:
                        verdict = "n/a (recombining in BCCH)"
                    else:
                        verdict = "model 1" if hi < 1 else "model 2" if lo > 1 else "unresolved"
                        if draws[zc]["_nblocks"] < 0:
                            verdict += " (window CI, too narrow)"
                    tests.append({"Z_compartment": zc, "autosome_set": ac, "stat": k, "ratio": pt,
                                  "lo": lo, "hi": hi, "favours": verdict})
        tests = pd.DataFrame(tests)
        tests.to_csv(f"{a.out}.test.tsv", sep="\t", index=False, float_format="%.4g")

    # console report
    pd.set_option("display.width", 200)
    show = summ[["compartment", "pair", "n_windows", "dxy", "dxy_n", "dxy_n_lo", "dxy_n_hi"] +
                (["years", "years_lo", "years_hi"] if "years" in summ else [])]
    print(show.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    if len(tests):
        print("\nZ BCCH-MOCH / autosomal BCCH-CACH (outgroup-scaled). Model 1: < 1; model 2: > 1")
        print(tests.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    plot(W, summ, comps, ingroup_pairs, f"{a.out}.png", out)
    print(f"\nwrote {a.out}.windows.tsv, .summary.tsv, .test.tsv, .png")


def plot(W, summ, comps, pairs, path, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#bbbbbb", "axes.labelcolor": "#191919",
                         "xtick.color": TEXT2, "ytick.color": TEXT2})
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(11, 8), gridspec_kw={"height_ratios": [1.3, 1]})

    # A: scaled dXY along the Z, with compartment bands
    shade = {"neo_PAR": "#f4f4f4", "MOCH_stratum": "#e2e2e2", "neo_nonPAR": "#c9c9c9",
             "anc_nonPAR": "#ffffff", "anc_nonPAR_CACHgap": "#ffffff", "uncertain": "#ffffff", "anc_PAR": "#f4f4f4"}
    for s, e, lab in Z_REGIONS:
        e = min(e, W[W.compartment.isin([r[2] for r in Z_REGIONS] + ["boundary"])].end.max())
        ax.axvspan(s / 1e6, e / 1e6, color=shade[lab], lw=0, zorder=0)
    for xm, lab in ((7.4, "neo-PAR"), (20.1, "MOCH stratum"), (43.8, "neo-PAR"), (64.1, "neo non-PAR"),
                    (108, "ancestral non-PAR")):
        ax.text(xm, 1.01, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8, color=TEXT2)
    zw = W[~W.compartment.str.startswith("autosomes")].sort_values("start")
    x = (zw.start + zw.end) / 2e6
    for p in pairs:
        y = zw[f"dxy_n:{p}"]
        ax.plot(x, y, color=PAIRCOL.get(p, "#888"), lw=0.5, alpha=0.35)
        seg = (zw.compartment != zw.compartment.shift()).cumsum()     # do not smooth across boundaries
        sm = y.groupby(seg).transform(lambda v: v.rolling(11, center=True, min_periods=3).median())
        ax.plot(x, sm, color=PAIRCOL.get(p, "#888"), lw=2, label=PAIRNAME.get(p, p))
        a = summ[(summ.compartment == "autosomes") & (summ.pair == p)]
        if len(a):
            ax.axhline(a.dxy_n.iloc[0], color=PAIRCOL.get(p, "#888"), lw=1, ls=(0, (4, 2)))
    ax.set_xlim(0, x.max())
    ax.set_xlabel("Position on Z (Mb, BCCH reference)")
    ax.set_ylabel(f"dXY / dXY to {out}")
    ax.set_title("A   Outgroup-scaled dXY along the Z (dashed = autosomal value)", loc="left", fontsize=12, pad=24)
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1, loc="lower left", fontsize=9)

    # B: compartment totals with block-bootstrap CIs
    k = len(pairs)
    for j, p in enumerate(pairs):
        s = summ[summ.pair == p].set_index("compartment").reindex(comps)
        xs = np.arange(len(comps)) + (j - (k - 1) / 2) * 0.22
        bx.errorbar(xs, s.dxy_n, yerr=[s.dxy_n - s.dxy_n_lo, s.dxy_n_hi - s.dxy_n], fmt="o", ms=6,
                    color=PAIRCOL.get(p, "#888"), ecolor=PAIRCOL.get(p, "#888"), elinewidth=1.5, capsize=0,
                    label=PAIRNAME.get(p, p))
    bx.set_xticks(range(len(comps)), [LABEL.get(c, c) for c in comps], rotation=20, ha="right")
    bx.set_ylabel(f"dXY / dXY to {out}")
    bx.set_title("B   Totals by compartment (95% block-bootstrap CI)", loc="left", fontsize=12)
    for a_ in (ax, bx):
        a_.grid(axis="y", color="#eeeeee", lw=0.8)
        a_.set_axisbelow(True)
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


if __name__ == "__main__":
    main()
