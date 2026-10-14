#!/usr/bin/env python3
# vcf2phylip_snps.py
# Converts biallelic SNP genotypes to a relaxed phylip alignment for raxml-ng
# heterozygotes are coded with IUPAC ambiguity codes, missing genotypes as N
# drops sites that are invariant across samples (required for raxml-ng ascertainment bias correction)
#
# usage: bcftools query -H -f '%REF\t%ALT[\t%GT]\n' in.vcf.gz | python3 vcf2phylip_snps.py out.phy

import sys

IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
         frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}

out_file = sys.argv[1]

header = sys.stdin.readline().rstrip("\n").split("\t")
# bcftools -H header looks like: # [1]REF  [2]ALT  [3]sample:GT ...
samples = [h.split("]", 1)[1].rsplit(":GT", 1)[0] for h in header[2:]]
seqs = [[] for _ in samples]

n_in = n_kept = n_invariant = n_skipped = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    n_in += 1
    ref, alt = f[0].upper(), f[1].upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        n_skipped += 1
        continue
    alleles = (ref, alt)
    het = IUPAC[frozenset((ref, alt))]
    col = []
    for gt in f[2:]:
        a = gt.replace("|", "/").split("/")
        if len(a) != 2 or "." in a:
            col.append("N")
        elif a[0] == a[1]:
            col.append(alleles[int(a[0])])
        else:
            col.append(het)
    # keep only columns with more than one state among called genotypes
    if len(set(c for c in col if c != "N")) < 2:
        n_invariant += 1
        continue
    for i, c in enumerate(col):
        seqs[i].append(c)
    n_kept += 1

width = max(len(s) for s in samples) + 2
with open(out_file, "w") as out:
    out.write(f"{len(samples)} {n_kept}\n")
    for s, seq in zip(samples, seqs):
        out.write(s.ljust(width) + "".join(seq) + "\n")

sys.stderr.write(f"Sites read: {n_in}\n"
                 f"Sites skipped (not biallelic ACGT SNPs): {n_skipped}\n"
                 f"Sites dropped (invariant across samples): {n_invariant}\n"
                 f"Sites written: {n_kept}\n")
