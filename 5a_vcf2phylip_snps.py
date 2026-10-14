#!/usr/bin/env python3
# vcf2phylip_snps.py
# Converts biallelic SNP genotypes to a relaxed phylip alignment for raxml-ng
# heterozygotes are coded with IUPAC ambiguity codes, missing genotypes as N
# drops sites that raxml-ng treats as invariant under ascertainment bias correction (ASC_LEWIS etc.):
#   raxml-ng expands ambiguity codes (R = A/G, ...) and calls a site invariant if one base is compatible
#   with every sample, e.g. A A R N is invariant. For a biallelic SNP that means a site is only kept if at
#   least one sample is homozygous for each allele.
#
# usage: bcftools query -H -f '%REF\t%ALT[\t%GT]\n' in.vcf.gz | python3 vcf2phylip_snps.py out.phy

import sys

STATES = {"A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
          "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"},
          "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"}}

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
    # raxml-ng invariant test: drop the column if any single base is compatible with every called genotype
    shared = {"A", "C", "G", "T"}
    for c in set(col):
        if c != "N":
            shared &= STATES[c]
    if shared:
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
                 f"Sites dropped (invariant for raxml-ng: no sample homozygous for one of the alleles): {n_invariant}\n"
                 f"Sites written: {n_kept}\n")
