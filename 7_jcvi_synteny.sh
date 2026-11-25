#!/bin/bash
# -----------------------------------------------------------------------------
# Whole-genome synteny (jcvi karyotype plots) between the BCCH, CACH and MOCH
# haplotype assemblies, and a table of every sequence syntenic to the BCCH sex
# chromosomes (W = CM059969.1, Z = CM059970.1; ancestral Z = 65,501,052-146,584,261,
# the rest of the Z = neo-Z).
#
# No gene annotation is needed: minimap2 assembly-to-assembly alignment blocks
# are used as synteny anchors (7a_paf2jcvi.py), and jcvi draws them at
# base-pair positions (--basepair).
#
# Alignments (minimap2 -x asm20, target = first genome of each pair):
#   bcch_h1 vs bcch_h2, cach_h1, cach_h2, moch_h1, moch_h2   (ordering + sex table)
#   cach_h1 vs cach_h2, moch_h1 vs moch_h2                   (haplotype comparisons)
# Plots (PDF + PNG in ${OUTDIR}/plots):
#   one per pair; species.h1 (cach_h1 / bcch_h1 / moch_h1);
#   all_haps (cach_h2 / cach_h1 / bcch_h1 / moch_h1 / moch_h2);
#   sex.<pair> (only sequences syntenic to the BCCH W/Z) for every bcch_h1 pair
# Ribbon colours: W = purple, ancestral Z = orange, neo-Z = green, other = grey.
#
# Usage:
#   micromamba activate jcvi
#   nohup bash 7_jcvi_synteny.sh > ${PROJDIR}/analyses/synteny/synteny.log 2>&1 &
#
# Optional env overrides:
#   THREADS=24  PRESET=asm20  MINLEN=1000000  MINBLOCK=20000  MINMAPQ=10
#   MIN_SEX_BP=500000  OUTDIR=...  SCRIPTDIR=...
#   Existing PAFs are reused; delete one to realign that pair.
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny
SCRIPTDIR=${SCRIPTDIR:-${PROJDIR}/scripts/Chickadee_NeoSex}
OUTDIR=${OUTDIR:-${PROJDIR}/analyses/synteny}
THREADS=${THREADS:-24}
PRESET=${PRESET:-asm20}
MINLEN=${MINLEN:-1000000}
MINBLOCK=${MINBLOCK:-20000}
MINMAPQ=${MINMAPQ:-10}
MIN_SEX_BP=${MIN_SEX_BP:-500000}
WCHR=CM059969.1
ZCHR=CM059970.1
ANC_START=65501052
ANC_END=146584261
REF_ACC=GCA_030490865.1

declare -A FASTA=(
  [bcch_h1]=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
  [bcch_h2]=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490855.1/GCA_030490855.1_bPoeAtr1.hap2_genomic.fna
  [cach_h1]=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap1.p_ctg.fa.gz
  [cach_h2]=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap2.p_ctg.fa.gz
  [moch_h1]=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/Poecile_gambeli_hap1_NCBI.fasta
  [moch_h2]=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/ncbi_hap2_submission/aspera_upload/Poecile_gambeli_hap2.fasta
)
GENOMES="bcch_h1 bcch_h2 cach_h1 cach_h2 moch_h1 moch_h2"
PAIRS="bcch_h1:bcch_h2 bcch_h1:cach_h1 bcch_h1:cach_h2 bcch_h1:moch_h1 bcch_h1:moch_h2 cach_h1:cach_h2 moch_h1:moch_h2"

mkdir -p ${OUTDIR}/sizes ${OUTDIR}/paf ${OUTDIR}/jcvi ${OUTDIR}/plots
cd ${OUTDIR}
date
echo "Script: 7_jcvi_synteny.sh   preset: ${PRESET}  threads: ${THREADS}"

for g in ${GENOMES}; do
  [[ -s ${FASTA[$g]} ]] || { echo "ERROR: ${g}: ${FASTA[$g]} missing"; exit 1; }
done

# --- 1. Sequence lengths (works for plain or gzipped FASTA) ---
for g in ${GENOMES}; do
  S=${OUTDIR}/sizes/${g}.sizes
  [[ -s ${S} ]] && continue
  zcat -f ${FASTA[$g]} | awk '/^>/ {if (n) print n"\t"l; n=substr($1,2); l=0; next} {l+=length($0)} END {if (n) print n"\t"l}' > ${S}
  echo "  ${g}: $(wc -l < ${S}) sequences, $(awk '$2>='${MINLEN}'' ${S} | wc -l) >= ${MINLEN} bp"
done

# --- 2. Chromosome names for BCCH (1, 1A, Z, W...) from the NCBI sequence report ---
NAMES=${OUTDIR}/sizes/bcch_h1.names.tsv
if [[ ! -s ${NAMES} ]] && command -v datasets > /dev/null && command -v dataformat > /dev/null; then
  datasets summary genome accession ${REF_ACC} --report sequence --as-json-lines \
    | dataformat tsv genome-seq --fields genbank-seq-acc,chr-name --elide-header > ${NAMES} \
    || { echo "  (could not fetch chromosome names; using accessions)"; rm -f ${NAMES}; }
fi

# --- 3. Whole-genome alignments ---
for pr in ${PAIRS}; do
  A=${pr%%:*}; B=${pr##*:}
  PAF=${OUTDIR}/paf/${A}.${B}.paf
  if [[ -s ${PAF} ]]; then echo "  ${A} vs ${B}: PAF exists"; continue; fi
  echo "[$(date)] minimap2 -x ${PRESET}: ${A} (target) vs ${B} (query)"
  minimap2 -x ${PRESET} --secondary=no -t ${THREADS} ${FASTA[$A]} ${FASTA[$B]} > ${PAF}.tmp
  mv ${PAF}.tmp ${PAF}
done

# --- 4. jcvi inputs + sex-synteny table ---
: > jcvi/genomes.tsv
for g in ${GENOMES}; do
  echo -e "${g}\t${OUTDIR}/sizes/${g}.sizes\t$([[ $g == bcch_h1 && -s ${NAMES} ]] && echo ${NAMES})" >> jcvi/genomes.tsv
done
: > jcvi/pairs.tsv
for pr in ${PAIRS}; do
  A=${pr%%:*}; B=${pr##*:}
  echo -e "${A}\t${B}\t${OUTDIR}/paf/${A}.${B}.paf" >> jcvi/pairs.tsv
done
{
  for pr in ${PAIRS}; do A=${pr%%:*}; B=${pr##*:}; echo -e "${A}.${B}\t${A},${B}\tall"; done
  for pr in ${PAIRS}; do A=${pr%%:*}; B=${pr##*:}; [[ $A == bcch_h1 ]] && echo -e "sex.${A}.${B}\t${A},${B}\tsex"; done
  echo -e "species.h1\tcach_h1,bcch_h1,moch_h1\tall"
  echo -e "sex.species.h1\tcach_h1,bcch_h1,moch_h1\tsex"
  echo -e "all_haps\tcach_h2,cach_h1,bcch_h1,moch_h1,moch_h2\tall"
} > jcvi/plots.tsv

python3 ${SCRIPTDIR}/7a_paf2jcvi.py --ref bcch_h1 \
  --genomes jcvi/genomes.tsv --pairs jcvi/pairs.tsv --plots jcvi/plots.tsv --outdir jcvi \
  --wchr ${WCHR} --zchr ${ZCHR} --anc-start ${ANC_START} --anc-end ${ANC_END} \
  --minlen ${MINLEN} --minblock ${MINBLOCK} --minmapq ${MINMAPQ} --min-sex-bp ${MIN_SEX_BP}

# --- 5. Karyotype plots ---
cd ${OUTDIR}/jcvi
cut -f1 plots.tsv | while read -r PLOT; do  # (stdin of jcvi redirected below)
  NT=$(wc -l < ${PLOT}.seqids)
  H=$(( NT < 3 ? 5 : 2 * NT ))
  for FMT in pdf png; do
    python -m jcvi.graphics.karyotype ${PLOT}.seqids ${PLOT}.layout --basepair --notex \
      --figsize 14x${H} --format ${FMT} --dpi 300 -o ${OUTDIR}/plots/${PLOT}.${FMT} < /dev/null \
      || echo "  WARNING: jcvi failed for ${PLOT}.${FMT}"
  done
done

echo
echo "Sequences flagged as syntenic to BCCH W/Z:"
awk -F'\t' 'NR==1 || $13=="yes"' ${OUTDIR}/jcvi/sex_syntenic_sequences.tsv | cut -f1-3,5-8,10,11
echo
echo "Plots: ${OUTDIR}/plots   Table: ${OUTDIR}/jcvi/sex_syntenic_sequences.tsv"
echo "Plot labels -> sequence names: ${OUTDIR}/jcvi/<genome>.alias.tsv"
date
