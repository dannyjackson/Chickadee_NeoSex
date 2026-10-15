#!/bin/bash
# B3_raxml.sh
# Maximum likelihood tree with raxml-ng from SNPs in one region of the original (unfiltered for PCA) Chr2 vcf
# region: CM059970.1:25440000-62170000
# filters: individuals by mean depth, genotypes by depth, sites by mean depth and missingness, optional thinning
# heterozygous genotypes are coded as IUPAC ambiguity codes; invariant sites are removed and
# ascertainment bias correction (Lewis) is applied, because the alignment contains SNPs only
# usage: sbatch B3_raxml.sh   (or bash B3_raxml.sh)
# run B1_PCA_setup.sh first to make the sample lists

#SBATCH --job-name=raxml_chr2
#SBATCH --output=slurm_output/raxml_chr2.%j
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=72:00:00

# activate the micromamba env with bcftools, raxml-ng, python, and R (ape, phangorn)
# install with: micromamba install -c conda-forge -c bioconda bcftools raxml-ng python r-base r-ape r-phangorn
ENV_NAME=your_env_name   # set to your env
eval "$(micromamba shell hook --shell bash)"
micromamba activate ${ENV_NAME}

if [ -z "${PROJDIR}" ]; then
    echo "Error: PROJDIR is not set" >&2
    exit 1
fi

CHR=Chr2
VCF=${PROJDIR}/data/vcfs/${CHR}.all_ind.diploid.vcf.gz
REGION=CM059970.1:25440000-62170000
REGIONNAME=CM059970.1_25.44-62.17Mb
SAMPLES=${PROJDIR}/reference_lists/all_samples.txt
SAMPLEINFO=${PROJDIR}/reference_lists/sample_info.txt
OUTDIR=${PROJDIR}/analyses/raxml/${CHR}/${REGIONNAME}
SCRIPTDIR=$(dirname "$(readlink -f "$0")")   # if submitting with sbatch, set this to the directory holding vcf2phylip_snps.py and B3_raxml_plot.r
THREADS=16
MIN_IND_DP=4       # drop individuals whose mean depth in the region is below this
MIN_GT_DP=3        # set genotypes with depth below this to missing
MIN_SITE_DP=4      # drop sites whose mean depth across kept individuals is below this
MAX_SITE_DP_X=2    # drop sites whose mean depth is above this multiple of the median individual depth (collapsed repeats/paralogs)
MAXMISS=0.2        # drop sites with more than this fraction of missing genotypes (after masking low-depth genotypes)
THIN=""            # optional: keep at most one SNP per this many bp (e.g. 1000) to shrink the alignment; empty = no thinning
BOOTSTRAPS=autoMRE{200}   # stop bootstrapping once support values converge (max 200); or set a number, e.g. 100
OUTGROUP=""        # optional, comma separated sample names to root the tree (e.g. BOCH samples); leave empty for midpoint rooting in the plot
PREFIX=${CHR}.${REGIONNAME}

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# region queries need an index on the input vcf
if [ ! -f ${VCF}.tbi ] && [ ! -f ${VCF}.csi ]; then
    echo "Error: no index for ${VCF}. Make one with: bcftools index -t ${VCF}" >&2
    exit 1
fi

# depth filters need per-genotype depth
if ! bcftools view -h ${VCF} | grep -q "^##FORMAT=<ID=DP,"; then
    echo "Error: ${VCF} has no FORMAT/DP field, so depth filters can't be applied" >&2
    exit 1
fi

# 1. mean depth per individual over biallelic snps in the region
bcftools view --threads ${THREADS} -r ${REGION} -S ${SAMPLES} --force-samples -m2 -M2 -v snps ${VCF} -Ou | \
bcftools stats -s - > ${PREFIX}.stats.txt

# PSC lines: column 3 = sample, column 10 = average depth (checked against the header below)
DPCOL=$(awk -F'\t' '$1 == "# PSC" {for (i = 1; i <= NF; i++) if ($i ~ /average depth/) print i; exit}' ${PREFIX}.stats.txt)
if [ "${DPCOL}" != "10" ]; then
    echo "Error: expected average depth in column 10 of bcftools stats PSC lines, found column '${DPCOL}'" >&2
    exit 1
fi
awk -F'\t' '$1 == "PSC" {print $3"\t"$10}' ${PREFIX}.stats.txt | sort -k2,2g > ${PREFIX}.ind_depth.txt

awk -v m=${MIN_IND_DP} '$2 >= m {print $1}' ${PREFIX}.ind_depth.txt > ${PREFIX}.samples_keep.txt
awk -v m=${MIN_IND_DP} '$2 <  m' ${PREFIX}.ind_depth.txt > ${PREFIX}.samples_dropped.txt

echo "Individuals kept (mean depth >= ${MIN_IND_DP}): $(wc -l < ${PREFIX}.samples_keep.txt)"
echo "Individuals dropped:"
cat ${PREFIX}.samples_dropped.txt

# median depth of kept individuals sets the upper site-depth cutoff
MEDIAN_DP=$(awk -v m=${MIN_IND_DP} '$2 >= m {print $2}' ${PREFIX}.ind_depth.txt | sort -g | \
    awk '{a[NR] = $1} END {if (NR % 2) print a[(NR + 1) / 2]; else print (a[NR / 2] + a[NR / 2 + 1]) / 2}')
MAX_SITE_DP=$(awk -v d=${MEDIAN_DP} -v x=${MAX_SITE_DP_X} 'BEGIN {print d * x}')
echo "Median individual depth: ${MEDIAN_DP}; site mean depth kept between ${MIN_SITE_DP} and ${MAX_SITE_DP}"

# 2. kept individuals only, biallelic snps
# 3. mask genotypes with depth < MIN_GT_DP
# 4. drop sites by mean depth, recompute AC/AN after masking, keep sites still variable, drop sites with too much missing data
bcftools view --threads ${THREADS} -r ${REGION} -S ${PREFIX}.samples_keep.txt ${VCF} -Ou | \
bcftools view --threads ${THREADS} -m2 -M2 -v snps -Ou | \
bcftools +setGT -Ou -- -t q -n . -i "FMT/DP < ${MIN_GT_DP}" | \
bcftools view --threads ${THREADS} -i "MEAN(FMT/DP) >= ${MIN_SITE_DP} && MEAN(FMT/DP) <= ${MAX_SITE_DP}" -Ou | \
bcftools +fill-tags -Ou -- -t AC,AN | \
bcftools view --threads ${THREADS} -i "AC > 0 && AC < AN && F_MISSING <= ${MAXMISS}" -Oz -o ${PREFIX}.snps.vcf.gz

bcftools index ${PREFIX}.snps.vcf.gz
echo "SNPs after depth and missingness filters: $(bcftools index -n ${PREFIX}.snps.vcf.gz)"

# optional thinning to one snp per THIN bp
if [ -n "${THIN}" ]; then
    bcftools +prune -w ${THIN}bp -n 1 -N 1st ${PREFIX}.snps.vcf.gz -Oz -o ${PREFIX}.snps.thin${THIN}.vcf.gz
    mv ${PREFIX}.snps.thin${THIN}.vcf.gz ${PREFIX}.snps.vcf.gz
    bcftools index -f ${PREFIX}.snps.vcf.gz
    echo "SNPs after thinning to 1 per ${THIN} bp: $(bcftools index -n ${PREFIX}.snps.vcf.gz)"
fi

# convert to phylip, heterozygotes as IUPAC codes, invariant columns removed
bcftools query -H -f '%REF\t%ALT[\t%GT]\n' ${PREFIX}.snps.vcf.gz | \
python3 ${SCRIPTDIR}/vcf2phylip_snps.py ${PREFIX}.phy

# check the alignment
raxml-ng --check --redo --msa ${PREFIX}.phy --model GTR+G+ASC_LEWIS --prefix ${PREFIX}.check

# compress to binary format for faster loading
raxml-ng --parse --redo --msa ${PREFIX}.phy --model GTR+G+ASC_LEWIS --prefix ${PREFIX}

# ml tree search + bootstraps, bootstrap support mapped onto the best tree
OUTGROUP_ARG=""
if [ -n "${OUTGROUP}" ]; then
    OUTGROUP_ARG="--outgroup ${OUTGROUP}"
fi

raxml-ng --all \
--msa ${PREFIX}.raxml.rba \
--model GTR+G+ASC_LEWIS \
--bs-trees ${BOOTSTRAPS} \
--threads ${THREADS} \
--seed 12345 \
${OUTGROUP_ARG} \
--prefix ${PREFIX}

# plot, tips colored by species; bootstrap tree: ${PREFIX}.raxml.support
Rscript ${SCRIPTDIR}/B3_raxml_plot.r ${OUTDIR}/${PREFIX}.raxml.support ${SAMPLEINFO} ${OUTDIR} ${PREFIX} "${OUTGROUP}"
