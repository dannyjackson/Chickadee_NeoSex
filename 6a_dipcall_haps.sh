#!/bin/bash
# -----------------------------------------------------------------------------
# Call phased diploid variants for ONE species directly from its two haplotype
# assemblies (hap1 + hap2 of the same individual) with dipcall.
#
# Both haplotypes are aligned to the shared BCCH reference (bPoeAtr1.hap1), so
# all species end up in the same coordinates for MSMC2 cross-coalescence.
#
# Outputs (in ${OUTBASE}/dipcall/<species>/):
#   <species>.dip.vcf.gz   phased VCF, GT = hap1|hap2 (FILTER "." = pass)
#   <species>.dip.bed      confident regions: exactly one >=50 kb alignment
#                          from EACH haplotype -> this is the callability mask
#   <species>.summary.txt  callable bp, het SNPs and heterozygosity per autosome
#
# Usage:
#   bash 6a_dipcall_haps.sh <species>      # bcch | cach | moch
#
# Run all three (from 6_MSMC_RCCR.md):
#   for sp in bcch cach moch; do
#     mkdir -p ${PROJDIR}/analyses/msmc_rccr/dipcall/${sp}
#     nohup bash ${SCRIPTDIR}/6a_dipcall_haps.sh ${sp} \
#       > ${PROJDIR}/analyses/msmc_rccr/dipcall/${sp}/dipcall.log 2>&1 &
#   done
#
# Optional env overrides:
#   REF=...  DIPCALL_DIR=~/programs/dipcall.kit  THREADS=16  MAKE_JOBS=2
#   ASM_PRESET=asm5   (minimap2 preset; dipcall default is asm5. If a divergent
#                      species, e.g. moch, has a much smaller dip.bed than bcch,
#                      rerun with ASM_PRESET=asm20 after deleting its outputs)
#   ZCHR=CM059970.1  WCHR=CM059969.1  MINLEN=1000000  (for the summary only)
# -----------------------------------------------------------------------------
set -euo pipefail

SP=${1:?Usage: $0 <species>   (bcch | cach | moch)}

PROJDIR=/nas/bcch_for_danny
REF=${REF:-${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna}
DIPCALL_DIR=${DIPCALL_DIR:-${HOME}/programs/dipcall.kit}
THREADS=${THREADS:-16}
MAKE_JOBS=${MAKE_JOBS:-2}          # hap1 and hap2 are aligned in parallel
ASM_PRESET=${ASM_PRESET:-asm5}
ZCHR=${ZCHR:-CM059970.1}           # Z (incl. neo-Z) - excluded
WCHR=${WCHR:-CM059969.1}           # W - excluded
MINLEN=${MINLEN:-1000000}          # drop unplaced scaffolds / mito
OUTBASE=${PROJDIR}/analyses/msmc_rccr
OUTDIR=${OUTBASE}/dipcall/${SP}

case ${SP} in
  bcch)
    HAP1=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
    HAP2=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490855.1/GCA_030490855.1_bPoeAtr1.hap2_genomic.fna ;;
  cach)
    HAP1=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap1.p_ctg.fa.gz
    HAP2=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap2.p_ctg.fa.gz ;;
  moch)
    HAP1=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/Poecile_gambeli_hap1_NCBI.fasta
    HAP2=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/ncbi_hap2_submission/aspera_upload/Poecile_gambeli_hap2.fasta ;;
  *) echo "ERROR: unknown species '${SP}' (expected bcch, cach or moch)"; exit 1 ;;
esac

mkdir -p ${OUTDIR}
cd ${OUTDIR}

echo "[$(date)] Species: ${SP}"
echo "  REF:  ${REF}"
echo "  HAP1: ${HAP1}"
echo "  HAP2: ${HAP2}"
echo "  minimap2 preset: ${ASM_PRESET}; threads: ${THREADS}"

for f in ${REF} ${HAP1} ${HAP2}; do
  [[ -s ${f} ]] || { echo "ERROR: ${f} missing or empty"; exit 1; }
done
[[ "${HAP1}" != "${HAP2}" ]] || { echo "ERROR: HAP1 and HAP2 are the same file"; exit 1; }
[[ -s ${REF}.fai ]] || { echo "ERROR: ${REF}.fai missing; run samtools faidx ${REF}"; exit 1; }

# --- 1. Build and run the dipcall makefile ---
# -a is REQUIRED: without it dipcall only keeps contigs named like human
# chromosomes (chr1..chr22/X), which would silently drop every CM0599xx.1 contig.
if [[ -s ${SP}.dip.vcf.gz && -s ${SP}.dip.bed ]]; then
  echo "[$(date)] dipcall outputs exist; skipping to post-processing"
else
  ${DIPCALL_DIR}/run-dipcall -a -t ${THREADS} ${SP} ${REF} ${HAP1} ${HAP2} > ${SP}.mak
  if [[ ${ASM_PRESET} != asm5 ]]; then
    sed -i "s/^MM2_OPT=-xasm5/MM2_OPT=-x${ASM_PRESET}/" ${SP}.mak
  fi
  echo "  minimap2 options used: $(grep '^MM2_OPT=' ${SP}.mak)"
  make -j${MAKE_JOBS} -f ${SP}.mak
fi

# --- 2. Name the sample after the species and index ---
if [[ ! -s ${SP}.dip.renamed.vcf.gz ]]; then
  echo ${SP} > sample_name.txt
  bcftools reheader -s sample_name.txt ${SP}.dip.vcf.gz -o ${SP}.dip.renamed.vcf.gz
  bcftools index -t -f ${SP}.dip.renamed.vcf.gz
fi

# --- 3. Sanity check: callable bp and heterozygosity per autosome ---
awk -v z=${ZCHR} -v w=${WCHR} -v m=${MINLEN} '$1!=z && $1!=w && $2>=m {print $1}' ${REF}.fai \
  > autosomes.txt
echo -e "chrom\tcallable_bp\tchrom_len\tfrac_callable\thet_snps\theterozygosity" > ${SP}.summary.txt
while read -r CHR; do
  LEN=$(awk -v c=${CHR} '$1==c {print $2}' ${REF}.fai)
  awk -v c=${CHR} -v OFS='\t' '$1==c' ${SP}.dip.bed > tmp.${CHR}.bed
  CALL=$(awk '{s+=$3-$2} END {print s+0}' tmp.${CHR}.bed)
  if [[ ${CALL} -gt 0 ]]; then
    HET=$(bcftools view -H -R tmp.${CHR}.bed -f .,PASS -v snps -m2 -M2 \
            -i 'GT="het"' ${SP}.dip.renamed.vcf.gz | wc -l)
  else
    HET=0
  fi
  awk -v c=${CHR} -v a=${CALL} -v l=${LEN} -v h=${HET} -v OFS='\t' \
    'BEGIN {print c, a, l, (l>0 ? a/l : 0), h, (a>0 ? h/a : 0)}' >> ${SP}.summary.txt
  rm -f tmp.${CHR}.bed
done < autosomes.txt

awk 'NR>1 {a+=$2; l+=$3; h+=$5} END {
  printf "[autosomes] callable %.0f of %.0f bp (%.1f%%); het SNPs %d; heterozygosity %.5f\n",
         a, l, 100*a/l, h, (a>0 ? h/a : 0)}' ${SP}.summary.txt

echo "[$(date)] Done with ${SP}"
