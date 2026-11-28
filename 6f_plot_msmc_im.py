#!/usr/bin/env python3
"""
6f_plot_msmc_im.py: summarise and plot MSMC-IM output (from 6c with RUN_IM=1, plus any extra beta runs).

    python3 ${SCRIPTDIR}/6f_plot_msmc_im.py ${PROJDIR}/analyses/msmc_rccr/output/msmcim \
        --out ${PROJDIR}/analyses/msmc_rccr/output/msmc_im

Reads every <pair>.b1_<x>.b2_<y>.MSMC_IM.estimates.txt in the directory. MSMC-IM reports time in
generations (it is given -mu); years = generations x the pair's mean generation length, as in 6e.

Per pair and per regularisation setting (beta1, beta2) it reports:
  * the times at which M(t) reaches 0.01, 0.25, 0.5, 0.75 and 0.99. M(t) is the cumulative probability that
    two lineages, one from each species, have merged into one population by time t (Wang et al. 2020).
    0.25-0.75 is the window Wang et al. use to summarise when separation happened.
  * the peak of the migration rate m(t) and when it occurs.
  * M(t) at --check-years (default 50 ka, 150 ka, 350 ka, 1 Ma) for a direct read of the BCCH-CACH plateau.

How to read it against the RCCR plateau (BCCH-CACH RCCR ~0.6-0.7 from ~0.35 to ~1.4 Ma):
  * Clean split: M(t) climbs from ~0 to ~1 over a short window and m(t) shows one hump near the split.
  * Old split plus later gene flow: going back in time, M(t) rises during the gene-flow period, stalls at an
    intermediate value (a shoulder) through the period of isolation, and reaches 1 only at the deep split.
    m(t) is raised during the gene-flow period, near zero through the isolation, and rises again at the
    deep split. The height of the shoulder approximates the share of ancestry exchanged.
  MSMC-IM fits symmetric, continuous migration, so a pulse is smeared over neighbouring intervals, and it
  cannot give the direction of gene flow. Compare beta settings: a feature that disappears when beta changes
  is not robust.

Outputs: <out>.summary.tsv, <out>.png
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

PAIRCOL = {"bcch_cach": "#6a4ce2", "bcch_moch": "#008a17", "cach_moch": "#ec4647"}   # as in 6e
NAMES = {"bcch": "Black-capped", "cach": "Carolina", "moch": "Mountain"}
GLACIAL = [("LGM", 21000, "#9d196a"), ("PGM", 140000, "#e78ac3")]
LEVELS = [0.01, 0.25, 0.5, 0.75, 0.99]


def parse_gen(s):
    return {k.lower(): float(v) for k, v in (x.split("=") for x in s.split(","))}


def read_est(fn):
    d = pd.read_csv(fn, sep=r"\s+")
    d = d.rename(columns={"left_boundaries": "left_time_boundary", "left_boundary": "left_time_boundary"})
    want = ["left_time_boundary", "im_N1", "im_N2", "m", "M"]
    if not set(want) <= set(d.columns):
        d = pd.read_csv(fn, sep=r"\s+", header=None, comment="#")
        if d.iloc[0].astype(str).str.contains("[A-Za-z]", regex=True).any():
            d = d.iloc[1:]
        d = d.iloc[:, :5].astype(float)
        d.columns = want
        print(f"note: {os.path.basename(fn)} has unexpected headers; assumed column order {want}", file=sys.stderr)
    return d[want].astype(float)


def cross_time(t, M, level):
    """First time (going back from the present) at which M(t) reaches `level`; log-time interpolation."""
    t, M = np.asarray(t), np.asarray(M)
    ok = t > 0
    t, M = t[ok], M[ok]
    idx = np.where(M >= level)[0]
    if not len(idx):
        return np.nan
    i = idx[0]
    if i == 0:
        return t[0]
    lt = np.interp(level, [M[i - 1], M[i]], np.log([t[i - 1], t[i]]))
    return float(np.exp(lt))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("imdir", help="directory with *.MSMC_IM.estimates.txt")
    ap.add_argument("--out", default="msmc_im")
    ap.add_argument("--gen", default="bcch=3.0,cach=2.3,moch=2.6", help="generation lengths (as in 6e)")
    ap.add_argument("--default-beta", default="1e-08,1e-06", help="beta1,beta2 drawn as the main line")
    ap.add_argument("--check-years", default="50000,150000,350000,1000000")
    a = ap.parse_args()
    gen = parse_gen(a.gen)
    checks = [float(x) for x in a.check_years.split(",")]
    b1d, b2d = (float(x) for x in a.default_beta.split(","))

    files = sorted(glob.glob(os.path.join(a.imdir, "*.MSMC_IM.estimates.txt")))
    if not files:
        sys.exit(f"no *.MSMC_IM.estimates.txt in {a.imdir}")
    runs, rows = [], []
    for f in files:
        m = re.match(r"(.+?)\.b1_([^.]+(?:\.[0-9]+)?(?:e-?\d+)?)\.b2_(.+?)\.MSMC_IM\.estimates\.txt$", os.path.basename(f))
        if not m:
            print(f"skip (name not recognised): {f}", file=sys.stderr)
            continue
        pair, b1, b2 = m.group(1), float(m.group(2)), float(m.group(3))
        sp = pair.split("_")
        g_mid = float(np.mean([gen[s] for s in sp]))
        g_lo, g_hi = min(gen[s] for s in sp), max(gen[s] for s in sp)
        d = read_est(f)
        if d.left_time_boundary.max() < 50:
            print(f"warning: {pair} times look coalescent-scaled, not generations; check MSMC-IM was run with -mu",
                  file=sys.stderr)
        d["years"] = d.left_time_boundary * g_mid
        r = {"pair": pair, "beta1": b1, "beta2": b2, "gen_mean": g_mid}
        for lv in LEVELS:
            tg = cross_time(d.left_time_boundary, d.M, lv)
            r[f"M{lv:g}_years"] = tg * g_mid
            r[f"M{lv:g}_years_lo"] = tg * g_lo
            r[f"M{lv:g}_years_hi"] = tg * g_hi
        pos = d[d.left_time_boundary > 0]
        # sanity check: a pair that is fully merged in the most recent intervals contradicts RCCR ~0 there,
        # which usually means the optimiser stalled (e.g. near its starting values)
        r["M_recent"] = float(pos.M.iloc[:3].mean())
        r["fit_flag"] = "CHECK FIT: M~1 at present" if r["M_recent"] > 0.5 else ""
        if r["fit_flag"]:
            print(f"warning: {pair} beta {b1:g},{b2:g}: M(t) = {r['M_recent']:.2f} in the most recent intervals; "
                  "the fit has probably failed. Check the fittingdetails PDF.", file=sys.stderr)
        k = pos.m.idxmax()
        r["peak_m"] = pos.m.max()
        r["peak_m_years"] = pos.loc[k, "years"]
        for c in checks:
            r[f"M_at_{c/1000:g}ka"] = float(np.interp(np.log(c), np.log(pos.years), pos.M))
        rows.append(r)
        runs.append((pair, b1, b2, d))

    S = pd.DataFrame(rows).sort_values(["pair", "beta1", "beta2"])
    S.to_csv(f"{a.out}.summary.tsv", sep="\t", index=False, float_format="%.6g")
    show = ["pair", "beta1", "beta2"] + [f"M{lv:g}_years" for lv in LEVELS] + ["peak_m_years"] + \
           [f"M_at_{c/1000:g}ka" for c in checks] + ["fit_flag"]
    pd.set_option("display.width", 250)
    print(S[show].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    plot(runs, b1d, b2d, a, f"{a.out}.png")
    print(f"\nwrote {a.out}.summary.tsv, {a.out}.png")


def plot(runs, b1d, b2d, a, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#bbbbbb", "axes.labelcolor": "#191919",
                         "xtick.color": "#666666", "ytick.color": "#666666"})
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    for pair, b1, b2, d in runs:
        d = d[d.left_time_boundary > 0]
        main = np.isclose(b1, b1d, rtol=1e-6, atol=0) and np.isclose(b2, b2d, rtol=1e-6, atol=0)
        col = PAIRCOL.get(pair, "#888888")
        lab = " vs ".join(NAMES.get(s, s) for s in pair.split("_")) if main else None
        kw = dict(color=col, lw=2 if main else 1, ls="-" if main else (0, (4, 2)), alpha=1 if main else 0.7,
                  where="post")
        ax.step(d.years, d.M, label=lab, **kw)
        bx.step(d.years, d.m, **kw)
    for x in (ax, bx):
        x.set_xscale("log")
        for name, yr, c in GLACIAL:
            x.axvline(yr, color=c, alpha=0.35, lw=7, zorder=0)
        x.grid(axis="y", color="#eeeeee", lw=0.8)
        x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    for lv in (0.25, 0.5, 0.75):
        ax.axhline(lv, color="#bbbbbb", lw=0.8, ls=":")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("M(t): cumulative merged probability")
    ax.set_title("A   MSMC-IM cumulative migration probability", loc="left", fontsize=12)
    bx.set_ylabel("m(t): symmetric migration rate\n(per generation)")
    bx.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    bx.set_title("B   Migration rate through time", loc="left", fontsize=12)
    bx.set_xlabel("Years before present (log scale; pair mean generation length)")
    h, l = ax.get_legend_handles_labels()
    h += [Line2D([], [], color="#666666", lw=2, label=f"beta {b1d:g}, {b2d:g}"),
          Line2D([], [], color="#666666", lw=1, ls=(0, (4, 2)), label="other beta settings")]
    h += [Line2D([], [], color=c, alpha=0.35, lw=7, label=n) for n, _, c in GLACIAL]
    ax.legend(handles=h, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(path, dpi=200)


if __name__ == "__main__":
    main()
