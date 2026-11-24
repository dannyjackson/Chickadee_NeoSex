#!/bin/bash
# -----------------------------------------------------------------------------
# Build MSMC2 multihetsep input files (one per autosome) jointly for all
# species, from the dipcall output of 6a_dipcall_haps.sh.
#
# For each species and autosome:
#   VCF  = PASS (FILTER ".") biallelic SNPs with a fully phased GT (hap1|hap2)
#   MASK = dipcall confident regions (dip.bed) for that chromosome
# then generate_multihetsep.py combines all species. Sites are only used where
# EVERY species' mask (plus the optional mappability / repeat masks) is callable.
#
# Haplotype indices in the output follow the order of SPECIES:
#   bcch = 0,1   cach = 2,3   moch = 4,5      (written to <tag>.indices.txt)
#
# Usage:
#   bash 6b_msmc_input.sh
#
# Optional env overrides:
#   SPECIES="bcch cach moch"     order sets the haplotype indices
#   REF=...  MSMCTOOLS=~/programs/msmc-tools
#   ZCHR=CM059970.1  WCHR=CM059969.1  MINLEN=1000000
#   AUTOSOMES=<file>             one chromosome name per line (default: all
#                                REF sequences >= MINLEN except Z and W)
#   MAPMASK_PATTERN=<path with {CHR}>   per-chromosome mappability mask, e.g.
#       ${PROJDIR}/analyses/snpable/bPoeAtr1_revised_{CHR}_mask.150.50.bed.gz
#   NEGMASK=<bed>                regions to EXCLUDE (e.g. RepeatMasker bed);
#                                split per chromosome automatically
# Outputs already present are reused; delete them to rebuild.
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny
REF=${REF:-${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna}
MSMCTOOLS=${MSMCTOOLS:-${HOME}/programs/msmc-tools}
SPECIES=${SPECIES:-"bcch cach moch"}
ZCHR=${ZCHR:-CM059970.1}
WCHR=${WCHR:-CM059969.1}
MINLEN=${MINLEN:-1000000}
MAPMASK_PATTERN=${MAPMASK_PATTERN:-}
NEGMASK=${NEGMASK:-}
OUTBASE=${PROJDIR}/analyses/msmc_rccr
TAG=$(echo ${SPECIES} | tr ' ' '_')

mkdir -p ${OUTBASE}/vcf ${OUTBASE}/mask ${OUTBASE}/input
cd ${OUTBASE}

# --- 0. Autosome list (Z carries the neo-Z region in this reference; drop Z and W) ---
AUTOSOMES=${AUTOSOMES:-${OUTBASE}/autosomes.txt}
if [[ ! -s ${AUTOSOMES} ]]; then
  awk -v z=${ZCHR} -v w=${WCHR} -v m=${MINLEN} '$1!=z && $1!=w && $2>=m {print $1}' ${REF}.fai \
    > ${AUTOSOMES}
fi
grep -qxF ${ZCHR} ${AUTOSOMES} && { echo "ERROR: Z (${ZCHR}) is in ${AUTOSOMES}"; exit 1; }
grep -qxF ${WCHR} ${AUTOSOMES} && { echo "ERROR: W (${WCHR}) is in ${AUTOSOMES}"; exit 1; }

echo "[$(date)] Species (index order): ${SPECIES}"
echo "  Autosomes: $(wc -l < ${AUTOSOMES}) from ${AUTOSOMES}"
echo "  Mappability mask: ${MAPMASK_PATTERN:-none}"
echo "  Negative mask:    ${NEGMASK:-none}"

# Record which haplotype indices belong to which species
i=0
: > ${OUTBASE}/input/${TAG}.indices.txt
for sp in ${SPECIES}; do
  echo -e "${sp}\t$((2*i)),$((2*i+1))" >> ${OUTBASE}/input/${TAG}.indices.txt
  i=$((i+1))
done
cat ${OUTBASE}/input/${TAG}.indices.txt

# --- 1. Per-species, per-chromosome VCFs and callability masks ---
for sp in ${SPECIES}; do
  DIPVCF=${OUTBASE}/dipcall/${sp}/${sp}.dip.renamed.vcf.gz
  DIPBED=${OUTBASE}/dipcall/${sp}/${sp}.dip.bed
  for f in ${DIPVCF} ${DIPBED}; do
    [[ -s ${f} ]] || { echo "ERROR: ${f} missing; run 6a_dipcall_haps.sh ${sp}"; exit 1; }
  done

  while read -r CHR; do
    VCF=${OUTBASE}/vcf/${sp}.${CHR}.vcf.gz
    MASK=${OUTBASE}/mask/${sp}.${CHR}.bed.gz
    if [[ ! -s ${VCF} ]]; then
      # GT regex is a safety net: PASS records from dipcall never have a '.' allele
      bcftools view -r ${CHR} -f .,PASS -v snps -m2 -M2 -Ov ${DIPVCF} \
        | awk 'BEGIN {OFS="\t"} /^#/ {print; next} $10 ~ /^[0-9]+\|[0-9]+(:|$)/' \
        | bgzip > ${VCF}
    fi
    if [[ ! -s ${MASK} ]]; then
      awk -v c=${CHR} -v OFS='\t' '$1==c {print $1, $2, $3}' ${DIPBED} | sort -k2,2n | bgzip > ${MASK}
    fi
  done < ${AUTOSOMES}
  echo "[$(date)] ${sp}: per-chromosome VCFs and masks ready"
done

# --- 2. Optional negative (e.g. repeat) mask, split per chromosome ---
if [[ -n ${NEGMASK} ]]; then
  mkdir -p ${OUTBASE}/mask/negative
  while read -r CHR; do
    NM=${OUTBASE}/mask/negative/neg.${CHR}.bed.gz
    [[ -s ${NM} ]] || awk -v c=${CHR} -v OFS='\t' '$1==c {print $1, $2, $3}' ${NEGMASK} \
                        | sort -k2,2n | bgzip > ${NM}
  done < ${AUTOSOMES}
fi

# --- 3. Joint multihetsep per chromosome ---
# Masks are read WITHOUT checking chromosome names, so every mask passed here
# must contain only this chromosome (done above).
while read -r CHR; do
  OUT=${OUTBASE}/input/msmc_input.${TAG}.${CHR}.txt
  if [[ -s ${OUT} ]]; then
    echo "  ${CHR}: input exists, skipping"; continue
  fi

  MASKS=""; VCFS=""
  for sp in ${SPECIES}; do
    MASKS="${MASKS} --mask=${OUTBASE}/mask/${sp}.${CHR}.bed.gz"
    VCFS="${VCFS} ${OUTBASE}/vcf/${sp}.${CHR}.vcf.gz"
  done
  if [[ -n ${MAPMASK_PATTERN} ]]; then
    MM=${MAPMASK_PATTERN//\{CHR\}/${CHR}}
    [[ -s ${MM} ]] || { echo "ERROR: mappability mask ${MM} missing"; exit 1; }
    MASKS="${MASKS} --mask=${MM}"
  fi
  [[ -n ${NEGMASK} ]] && MASKS="${MASKS} --negative_mask=${OUTBASE}/mask/negative/neg.${CHR}.bed.gz"

  python3 ${MSMCTOOLS}/generate_multihetsep.py ${MASKS} ${VCFS} > ${OUT}

  # drop empty files so they never reach msmc2
  if [[ ! -s ${OUT} ]]; then
    echo "  WARNING: ${CHR} produced no segregating sites; removed"; rm -f ${OUT}
  else
    echo "  ${CHR}: $(wc -l < ${OUT}) segregating sites"
  fi
done < ${AUTOSOMES}

ls ${OUTBASE}/input/msmc_input.${TAG}.*.txt > ${OUTBASE}/input/${TAG}.input_files.txt
echo "[$(date)] Wrote $(wc -l < ${OUTBASE}/input/${TAG}.input_files.txt) input files; list in input/${TAG}.input_files.txt"
