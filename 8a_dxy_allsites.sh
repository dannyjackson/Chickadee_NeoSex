#!/bin/bash
# 8a_dxy_allsites.sh
# All-sites calls (males only) on the Z and autosomes, then windowed pi and dXY with pixy.
#
# Why a new call: the existing VCFs (1a_Call_Variants.md) were made with `bcftools call -mv`,
# i.e. variant sites only. dXY needs every callable site in the denominator, otherwise it is
# inflated, and by different amounts on the Z and the autosomes.
#
# Males only: males are ZZ, so the Z is called as a true diploid; females would appear as
# false homozygotes on the Z. The same males are used for the autosomes so the compartments
# are directly comparable. BOCH males are the outgroup used in 8b to scale mutation rate.
#
# Usage (all variables can be overridden on the command line):
#   nohup bash ${SCRIPTDIR}/8a_dxy_allsites.sh > ${PROJDIR}/analyses/dxy/8a.log 2>&1 &
#   DO_AUTOS=0 bash 8a_dxy_allsites.sh                    # Z only (run this first; much faster)
#   AUTOS=my_autosome_subset.txt bash 8a_dxy_allsites.sh  # a subset of autosomes, one name per line
#   NEGMASK=repeats.bed bash 8a_dxy_allsites.sh           # exclude repeat / low-mappability regions
#
# Re-running skips chunks that already finished, so an interrupted run can be restarted.

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate ${ENV:-bcf}            # needs bcftools, htslib and pixy (see 8_dXY.md)

PROJDIR=${PROJDIR:-/nas/bcch_for_danny}
REF_noW=${REF_noW:-${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.no_W.fna}
BAMLIST=${BAMLIST:-${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt}
SEXDIR=${SEXDIR:-${PROJDIR}/analyses/sex_inference}
AUTOS=${AUTOS:-${PROJDIR}/reference_lists/bcch_ref.autosomes.txt}
OUTDIR=${OUTDIR:-${PROJDIR}/analyses/dxy}

ZCHR=${ZCHR:-CM059970.1}
SPECIES=${SPECIES:-"bcch cach moch boch"}
DO_AUTOS=${DO_AUTOS:-1}
THREADS=${THREADS:-24}           # shared server: max 24 cores
CHUNK=${CHUNK:-5000000}          # calling chunk size (bp); one chunk per core at a time
WIN=${WIN:-100000}               # pixy window size; 100 kb matches the Twisst windows
MAPQ=${MAPQ:-20}
BASEQ=${BASEQ:-20}
MIN_DP=${MIN_DP:-5}              # genotypes with DP below this are set to missing
MAX_MISS=${MAX_MISS:-0.5}        # drop sites missing in more than this fraction of samples
NEGMASK=${NEGMASK:-}             # optional BED of regions to exclude

mkdir -p ${OUTDIR}/{lists,chunks,vcf,pixy}
cd ${OUTDIR}

echo "[$(date)] 8a_dxy_allsites.sh"
echo "  Reference:   ${REF_noW}"
echo "  Z:           ${ZCHR}; autosomes: $([[ ${DO_AUTOS} == 1 ]] && echo ${AUTOS} || echo skipped)"
echo "  Threads:     ${THREADS}; chunk ${CHUNK}; window ${WIN}"
echo "  Filters:     MAPQ>=${MAPQ}, BASEQ>=${BASEQ}, genotype DP>=${MIN_DP}, site missing<=${MAX_MISS}"
echo "  Neg. mask:   ${NEGMASK:-none}"

## 1. Male sample lists and the pixy populations file
: > lists/males.txt
: > lists/populations.txt
for sp in ${SPECIES}; do
    SEX=${SEXDIR}/${sp}/sex_calls.tsv
    [[ -s ${SEX} ]] || { echo "Missing ${SEX}; run 2_Infer_Sex.md first" >&2; exit 1; }
    SPUP=$(echo ${sp} | tr '[:lower:]' '[:upper:]')
    # sex_calls.tsv: header line, sample in column 1, sex (M/F) in the last column
    awk 'NR > 1 && $NF == "M" {print $1}' ${SEX} > lists/${sp}.males.txt
    awk -v sp=${SPUP} '{print $1"\t"sp}' lists/${sp}.males.txt >> lists/populations.txt
    cat lists/${sp}.males.txt >> lists/males.txt
    echo "  ${SPUP}: $(wc -l < lists/${sp}.males.txt) males"
done

## 2. Region chunks (Z first, then autosomes)
[[ -s ${REF_noW}.fai ]] || samtools faidx ${REF_noW}
{ echo ${ZCHR}; [[ ${DO_AUTOS} == 1 ]] && cat ${AUTOS}; } > lists/chroms.txt
awk -v c=${CHUNK} 'NR == FNR {keep[$1] = 1; next}
    ($1 in keep) {for (s = 1; s <= $2; s += c) {e = s + c - 1; if (e > $2) e = $2; print $1":"s"-"e}}' \
    lists/chroms.txt ${REF_noW}.fai > lists/regions.txt
echo "  $(wc -l < lists/regions.txt) chunks across $(wc -l < lists/chroms.txt) chromosomes"

## 3. All-sites calls per chunk, in parallel
call_chunk() {
    set -o pipefail
    local region=$1
    local out=chunks/$(echo ${region} | tr ':-' '__').vcf.gz
    [[ -s ${out}.tbi ]] && return 0                       # finished in an earlier run
    local mask=()
    [[ -n ${NEGMASK} ]] && mask=(-T ^${NEGMASK})
    bcftools mpileup -f ${REF_noW} -r ${region} -b ${BAMLIST} -S lists/males.txt \
        -q ${MAPQ} -Q ${BASEQ} -a AD,DP -Ou |
    bcftools call -m -Ou |                                  # -m without -v: keep invariant sites
    bcftools view -M2 -V indels,mnps,other ${mask[@]+"${mask[@]}"} -Ou |   # invariant + biallelic SNPs only
    bcftools filter -S . -e "FMT/DP < ${MIN_DP}" -Ou |      # low-depth genotypes -> missing
    bcftools view -e "TYPE=\"snp\" && QUAL < 20" -Ou |
    bcftools view -i "F_MISSING <= ${MAX_MISS}" -Oz -o ${out}.tmp.vcf.gz
    mv ${out}.tmp.vcf.gz ${out}
    tabix -f -p vcf ${out}
    echo "[$(date)] done ${region}"
}
export -f call_chunk
export REF_noW BAMLIST MAPQ BASEQ MIN_DP MAX_MISS NEGMASK

echo "[$(date)] Calling (${THREADS} chunks at a time)"
xargs -a lists/regions.txt -P ${THREADS} -I {} bash -c 'call_chunk "$@"' _ {}

## 4. Concatenate: one VCF for the Z, one for the autosomes
concat() {      # $1 = output, remaining = chromosomes
    local out=$1; shift
    [[ -s ${out}.tbi ]] && { echo "  exists: ${out}"; return 0; }
    for c in "$@"; do grep "^${c}:" lists/regions.txt; done |
        awk '{gsub(/[:-]/, "_"); print "chunks/"$0".vcf.gz"}' > ${out}.files
    bcftools concat -f ${out}.files --threads ${THREADS} -Oz -o ${out}
    tabix -f -p vcf ${out}                                  # pixy needs a .tbi index
}
concat vcf/allsites.Z.males.vcf.gz ${ZCHR}
[[ ${DO_AUTOS} == 1 ]] && concat vcf/allsites.autosomes.males.vcf.gz $(cat ${AUTOS})

## 5. pixy: windowed pi (per species) and dXY (all species pairs)
run_pixy() {    # $1 = vcf, $2 = prefix
    [[ -s pixy/${2}_dxy.txt ]] && { echo "  exists: pixy/${2}_dxy.txt"; return 0; }
    pixy --stats pi dxy --vcf $1 --populations lists/populations.txt \
         --window_size ${WIN} --n_cores ${THREADS} \
         --output_folder pixy --output_prefix $2
}
echo "[$(date)] pixy"
run_pixy vcf/allsites.Z.males.vcf.gz Z
[[ ${DO_AUTOS} == 1 ]] && run_pixy vcf/allsites.autosomes.males.vcf.gz autosomes

echo "[$(date)] Done. Next: python3 8b_dxy_summary.py ${OUTDIR}/pixy"
