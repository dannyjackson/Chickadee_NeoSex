#!/usr/bin/env python3
"""
Convert minimap2 whole-genome PAFs into jcvi karyotype inputs, and list every
sequence syntenic to the reference sex chromosomes.

Each PAF alignment block (primary, MAPQ >= --minmapq, length >= --minblock)
becomes one synteny ribbon: two point features (block start / end) per genome in
<genome>.bed and one line in <A>.<B>.simple. jcvi is run with --basepair, so
ribbons are drawn at true base-pair positions.

Ribbons are coloured by what the reference-side (A) interval overlaps:
  W            -> --wcolor      (whole W)
  ancestral Z  -> --anczcolor   (Z between --anc-start and --anc-end)
  neo-Z        -> --neozcolor   (rest of Z)
For pairs that do not include the reference (e.g. cach_h1 vs cach_h2), the
sex-linked intervals of A are taken from A's own alignment to the reference.

Inputs
  --genomes TSV : label <tab> sizes file (name<tab>length, file order) [<tab> alias map]
  --pairs   TSV : A <tab> B <tab> PAF   (PAF made with: minimap2 A.fa B.fa, so A = target)
  --plots   TSV : plot_name <tab> track1,track2,... (top to bottom) <tab> all|sex
Outputs (in --outdir)
  <genome>.bed, <A>.<B>.simple, <plot>.seqids, <plot>.layout,
  <genome>.alias.tsv            alias used on the plot <-> real sequence name
  sex_syntenic_sequences.tsv    every sequence with alignment to ref W / Z
"""
import argparse
import collections
import os
import sys

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--ref", required=True, help="label of the reference genome (BCCH hap1)")
p.add_argument("--genomes", required=True)
p.add_argument("--pairs", required=True)
p.add_argument("--plots", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--wchr", default="CM059969.1")
p.add_argument("--zchr", default="CM059970.1")
p.add_argument("--anc-start", type=int, default=65501052)
p.add_argument("--anc-end", type=int, default=146584261)
p.add_argument("--minlen", type=int, default=1000000, help="min sequence length to plot")
p.add_argument("--minblock", type=int, default=20000, help="min alignment block (bp) to draw")
p.add_argument("--minblock-table", type=int, default=5000, help="min block (bp) counted in the sex table")
p.add_argument("--minmapq", type=int, default=10)
p.add_argument("--min-plot-aln", type=int, default=100000,
               help="a sequence is plotted only if >= this many bp align to plotted ref chromosomes")
p.add_argument("--min-sex-bp", type=int, default=500000, help="bp aligned to W/Z to flag a sequence")
p.add_argument("--wcolor", default="#7570b3")
p.add_argument("--anczcolor", default="#d95f02")
p.add_argument("--neozcolor", default="#1b9e77")
a = p.parse_args()

os.makedirs(a.outdir, exist_ok=True)
CLASSES = ("W", "ancZ", "neoZ")
COLOR = {"W": a.wcolor, "ancZ": a.anczcolor, "neoZ": a.neozcolor}


def read_tsv(fn):
    with open(fn) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]


# ---------- genomes ----------
sizes, order_in_file, alias_in = {}, {}, {}
for row in read_tsv(a.genomes):
    g, sz = row[0], row[1]
    sizes[g], order_in_file[g] = {}, []
    for name, ln in read_tsv(sz):
        sizes[g][name] = int(ln)
        order_in_file[g].append(name)
    alias_in[g] = {}
    if len(row) > 2 and row[2] and os.path.exists(row[2]):
        for r in read_tsv(row[2]):
            if len(r) >= 2 and r[1]:
                alias_in[g][r[0]] = r[1]
if a.ref not in sizes:
    sys.exit(f"--ref {a.ref} not in {a.genomes}")


def ref_overlap(chrom, s, e):
    out = {}
    if chrom == a.wchr:
        out["W"] = e - s
    elif chrom == a.zchr:
        anc = max(0, min(e, a.anc_end) - max(s, a.anc_start - 1))
        if anc:
            out["ancZ"] = anc
        if e - s - anc > 0:
            out["neoZ"] = e - s - anc
    return out


# ---------- read PAFs ----------
def read_paf(fn, minblock):
    blocks = []
    with open(fn) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = f[12:]
            if "tp:A:S" in tags or "tp:A:i" in tags:
                continue
            if int(f[11]) < a.minmapq or int(f[10]) < minblock:
                continue
            # query = B, target = A
            blocks.append(dict(b=f[0], bs=int(f[2]), be=int(f[3]), strand=f[4],
                               a=f[5], as_=int(f[7]), ae=int(f[8]), alen=int(f[10])))
    return blocks


pairs = read_tsv(a.pairs)
pair_blocks = {}
for A, B, paf in pairs:
    if not os.path.exists(paf):
        sys.exit(f"missing PAF {paf}")
    pair_blocks[(A, B)] = read_paf(paf, min(a.minblock, a.minblock_table))
    print(f"[paf2jcvi] {A} vs {B}: {len(pair_blocks[(A, B)])} blocks", file=sys.stderr)

# ---------- sex-linked intervals and table (from alignments to the reference) ----------
sex_iv = collections.defaultdict(list)  # genome -> [(chrom, s, e, class)]
for c in (a.wchr, a.zchr):
    if c in sizes[a.ref]:
        L = sizes[a.ref][c]
        if c == a.wchr:
            sex_iv[a.ref].append((c, 0, L, "W"))
        else:
            sex_iv[a.ref].append((c, 0, a.anc_start - 1, "neoZ"))
            sex_iv[a.ref].append((c, a.anc_start - 1, a.anc_end, "ancZ"))
            sex_iv[a.ref].append((c, a.anc_end, L, "neoZ"))

table_rows = []
best_ref = {}  # genome -> {chrom: (ref chrom, median ref pos, minus fraction, aligned bp to plotted ref)}
ref_plot = [c for c in order_in_file[a.ref] if sizes[a.ref][c] >= a.minlen]
ref_rank = {c: i for i, c in enumerate(ref_plot)}

for (A, B), blocks in pair_blocks.items():
    if A != a.ref:
        continue
    per = collections.defaultdict(lambda: dict(total=0, W=0, ancZ=0, neoZ=0, auto=collections.Counter(),
                                               zmin=None, zmax=None, wmin=None, wmax=None))
    hits = collections.defaultdict(lambda: collections.defaultdict(list))  # b -> refchrom -> [(mid, alen, strand)]
    for k in blocks:
        if k["alen"] >= a.minblock_table:
            d = per[k["b"]]
            d["total"] += k["alen"]
            ov = ref_overlap(k["a"], k["as_"], k["ae"])
            if ov:
                span = k["ae"] - k["as_"]
                for cl, bp in ov.items():
                    d[cl] += round(k["alen"] * bp / span)  # split by ref overlap (approximate)
                key = "z" if k["a"] == a.zchr else "w"
                d[key + "min"] = k["as_"] if d[key + "min"] is None else min(d[key + "min"], k["as_"])
                d[key + "max"] = k["ae"] if d[key + "max"] is None else max(d[key + "max"], k["ae"])
                main = max(ov, key=ov.get)
                sex_iv[B].append((k["b"], k["bs"], k["be"], main))
            else:
                d["auto"][k["a"]] += k["alen"]
        if k["alen"] >= a.minblock and k["a"] in ref_rank:
            hits[k["b"]][k["a"]].append(((k["as_"] + k["ae"]) / 2, k["alen"], k["strand"]))

    best_ref[B] = {}
    for b, byref in hits.items():
        tot = {r: sum(x[1] for x in v) for r, v in byref.items()}
        r = max(tot, key=tot.get)
        v = sorted(byref[r])
        half, acc, med = tot[r] / 2, 0, v[0][0]
        for mid, ln, _ in v:
            acc += ln
            if acc >= half:
                med = mid
                break
        minus = sum(ln for _, ln, st in byref[r] if st == "-") / tot[r]
        best_ref[B][b] = (r, med, minus, sum(tot.values()))

    for b, d in per.items():
        sexbp = d["W"] + d["ancZ"] + d["neoZ"]
        if sexbp == 0:
            continue
        L = sizes[B].get(b, 0)
        top_auto = d["auto"].most_common(1)
        flag = sexbp >= a.min_sex_bp or (sexbp >= 100000 and sexbp / max(d["total"], 1) >= 0.5)
        table_rows.append([B, b, L, d["total"], d["W"], d["ancZ"], d["neoZ"], sexbp,
                           f"{sexbp / max(d['total'], 1):.3f}",
                           f"{d['zmin']}-{d['zmax']}" if d["zmin"] is not None else "NA",
                           f"{d['wmin']}-{d['wmax']}" if d["wmin"] is not None else "NA",
                           f"{top_auto[0][0]}:{top_auto[0][1]}" if top_auto else "NA",
                           "yes" if flag else "no"])

table_rows.sort(key=lambda r: (r[0], -r[7]))
with open(os.path.join(a.outdir, "sex_syntenic_sequences.tsv"), "w") as fh:
    fh.write("\t".join(["genome", "sequence", "length", "aligned_to_ref_bp", "bp_to_W", "bp_to_ancestral_Z",
                        "bp_to_neo_Z", "bp_to_sex_total", "frac_aligned_that_is_sex",
                        "span_on_ref_Z", "span_on_ref_W", "top_autosome_hit:bp", "flagged"]) + "\n")
    for r in table_rows:
        fh.write("\t".join(map(str, r)) + "\n")

# ---------- which sequences to plot, order, orientation, alias ----------
plot_seqs, rev, alias = {}, {}, {}
plot_seqs[a.ref], rev[a.ref] = ref_plot, set()
for g in sizes:
    if g == a.ref:
        continue
    br = best_ref.get(g, {})
    keep = [b for b, (r, med, minus, tot) in br.items()
            if sizes[g].get(b, 0) >= a.minlen and tot >= a.min_plot_aln]
    keep.sort(key=lambda b: (ref_rank[br[b][0]], br[b][1]))
    plot_seqs[g] = keep
    rev[g] = {b for b in keep if br[b][2] > 0.5}
    if not br:
        print(f"[paf2jcvi] WARNING: {g} has no alignment to {a.ref}; it will not be plotted", file=sys.stderr)


def mk_alias(g):
    out, used = {}, collections.Counter()
    if g == a.ref:
        for c in plot_seqs[g]:
            nm = alias_in[g].get(c) or ("W" if c == a.wchr else "Z" if c == a.zchr else c)
            out[c] = nm
    else:
        refal = alias[a.ref]
        for b in plot_seqs[g]:
            base = refal[best_ref[g][b][0]]
            used[base] += 1
            out[b] = base if used[base] == 1 else f"{base}{chr(ord('a') + used[base] - 1)}"
    # aliases must be unique and comma/space free
    clean = {}
    seen = set()
    for k, v in out.items():
        v = v.replace(",", "_").replace(" ", "_")
        while v in seen:
            v += "_"
        seen.add(v)
        clean[k] = v
    return clean


alias[a.ref] = mk_alias(a.ref)
for g in sizes:
    if g != a.ref:
        alias[g] = mk_alias(g)
for g in sizes:
    with open(os.path.join(a.outdir, f"{g}.alias.tsv"), "w") as fh:
        fh.write("alias\tsequence\tlength\treversed_on_plot\tbest_ref_chrom\n")
        for s in plot_seqs[g]:
            br = best_ref.get(g, {}).get(s)
            fh.write(f"{alias[g][s]}\t{s}\t{sizes[g][s]}\t{'yes' if s in rev[g] else 'no'}\t"
                     f"{br[0] if br else s}\n")

# ---------- beds + simple files ----------
beds = collections.defaultdict(list)
for g in sizes:
    for s in plot_seqs[g]:  # endpoint features so jcvi draws full-length chromosomes
        beds[g].append((alias[g][s], 0, 1, f"{g}|{s}|start"))
        beds[g].append((alias[g][s], sizes[g][s] - 1, sizes[g][s], f"{g}|{s}|end"))


sex_by_chrom = collections.defaultdict(list)
for g_, ivs in sex_iv.items():
    for c, xs, xe, cl in ivs:
        sex_by_chrom[(g_, c)].append((xs, xe, cl))


def sex_class(g, chrom, s, e):
    """class covering >= 50% of the interval (overlaps summed per class)"""
    bp = collections.Counter()
    for xs, xe, cl in sex_by_chrom.get((g, chrom), []):
        ov = min(e, xe) - max(s, xs)
        if ov > 0:
            bp[cl] += ov
    if not bp:
        return None
    cl, v = bp.most_common(1)[0]
    return cl if v >= 0.5 * (e - s) else None


for (A, B), blocks in pair_blocks.items():
    pa, pb = set(plot_seqs[A]), set(plot_seqs[B])
    lines, n = [], 0
    for k in blocks:
        if k["alen"] < a.minblock or k["a"] not in pa or k["b"] not in pb:
            continue
        n += 1
        ida, idb = f"{A}.{B}.{n}", f"{B}.{A}.{n}"
        beds[A] += [(alias[A][k["a"]], k["as_"], k["as_"] + 1, ida + "s"),
                    (alias[A][k["a"]], k["ae"] - 1, k["ae"], ida + "e")]
        beds[B] += [(alias[B][k["b"]], k["bs"], k["bs"] + 1, idb + "s"),
                    (alias[B][k["b"]], k["be"] - 1, k["be"], idb + "e")]
        cl = sex_class(A, k["a"], k["as_"], k["ae"])
        hl = f"{COLOR[cl]}*" if cl else ""
        lines.append(f"{hl}{ida}s\t{ida}e\t{idb}s\t{idb}e\t{k['alen']}\t{k['strand']}")
    with open(os.path.join(a.outdir, f"{A}.{B}.simple"), "w") as fh:
        fh.write("\n".join(lines) + ("\n" if lines else ""))
    print(f"[paf2jcvi] {A}.{B}.simple: {len(lines)} ribbons", file=sys.stderr)

for g, rows in beds.items():
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(os.path.join(a.outdir, f"{g}.bed"), "w") as fh:
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")

# ---------- seqids + layout per plot ----------
pairset = {(A, B) for A, B, _ in pairs}
sex_seqs = {g: {c for c, _, _, _ in sex_iv.get(g, [])} for g in sizes}
for name, tracks, mode in read_tsv(a.plots):
    tr = tracks.split(",")
    n = len(tr)
    ys = [0.85 - i * (0.7 / max(n - 1, 1)) for i in range(n)] if n > 1 else [0.5]
    seqid_lines = []
    for g in tr:
        ss = plot_seqs[g]
        if mode == "sex":
            ss = [s for s in ss if s in sex_seqs[g]]
        if not ss:
            sys.exit(f"plot {name}: no sequences to draw for {g} (mode {mode})")
        seqid_lines.append(",".join(alias[g][s] + ("-" if s in rev[g] else "") for s in ss))
    with open(os.path.join(a.outdir, f"{name}.seqids"), "w") as fh:
        fh.write("\n".join(seqid_lines) + "\n")
    with open(os.path.join(a.outdir, f"{name}.layout"), "w") as fh:
        fh.write("# y, xstart, xend, rotation, color, label, va, bed\n")
        for i, g in enumerate(tr):
            va = "top" if i < n / 2 else "bottom"
            fh.write(f"{ys[i]:.3f}, .1, .9, 0, , {g}, {va}, {g}.bed\n")
        fh.write("# edges\n")
        for i in range(n - 1):
            g1, g2 = tr[i], tr[i + 1]
            if (g1, g2) in pairset:
                fh.write(f"e, {i}, {i + 1}, {g1}.{g2}.simple\n")
            elif (g2, g1) in pairset:
                fh.write(f"e, {i + 1}, {i}, {g2}.{g1}.simple\n")
            else:
                print(f"[paf2jcvi] WARNING: plot {name}: no alignment between {g1} and {g2}", file=sys.stderr)

flagged = [r for r in table_rows if r[-1] == "yes"]
print(f"[paf2jcvi] {len(flagged)} sequences flagged as syntenic to ref W/Z "
      f"(see sex_syntenic_sequences.tsv)", file=sys.stderr)
for r in flagged:
    print(f"  {r[0]}\t{r[1]}\tW={r[4]}\tancZ={r[5]}\tneoZ={r[6]}", file=sys.stderr)
