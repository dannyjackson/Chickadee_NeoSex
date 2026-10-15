# B3_raxml.sh
# Maximum likelihood tree with raxml-ng from SNPs in one region of the original (unfiltered for PCA) Chr2 vcf
# region: CM059970.1:25440000-62170000
# heterozygous genotypes are coded as IUPAC ambiguity codes; invariant sites are removed and
# ascertainment bias correction (Lewis) is applied, because the alignment contains SNPs only
# usage: sbatch B3_raxml.sh   (or bash B3_raxml.sh)
# run B1_PCA_setup.sh first to make the sample lists


# activate the micromamba env with bcftools, raxml-ng, python, and R (ape, phangorn)
# install with: micromamba install -c conda-forge -c bioconda bcftools raxml-ng python r-base r-ape r-phangorn

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd ${SCRIPTDIR}
git pull
REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna

CHR=NeoZ
VCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz 
OUTDIR=${PROJDIR}/analyses/raxml/${CHR}


CHR=NeoZ
VCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz
REGION=CM059970.1:25440000-62170000
REGIONNAME=CM059970.1_25.44-62.17Mb
SAMPLES=${PROJDIR}/reference_lists/all_samples.txt
SAMPLEINFO=${PROJDIR}/reference_lists/sample_info.txt
OUTDIR=${PROJDIR}/analyses/raxml/${CHR}/${REGIONNAME}
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
fi

# 1. mean depth per individual over biallelic snps in the region
bcftools view --threads ${THREADS} -r ${REGION} -S ${SAMPLES} --force-samples -m2 -M2 -v snps ${VCF} -Ou | \
bcftools stats -s - > ${PREFIX}.stats.txt

# PSC lines: column 3 = sample, column 10 = average depth (checked against the header below)
DPCOL=$(awk -F'\t' '$1 == "# PSC" {for (i = 1; i <= NF; i++) if ($i ~ /average depth/) print i; exit}' ${PREFIX}.stats.txt)
if [ "${DPCOL}" != "10" ]; then
    echo "Error: expected average depth in column 10 of bcftools stats PSC lines, found column '${DPCOL}'" >&2
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

# subset region and samples, keep biallelic snps that are variable in these samples, drop sites with too much missing data
bcftools view --threads ${THREADS} -r ${REGION} -S ${PREFIX}.samples_keep.txt ${VCF} -Ou | \
bcftools view --threads ${THREADS} -m2 -M2 -v snps -Ou | \
bcftools +setGT -Ou -- -t q -n . -i "FMT/DP < ${MIN_GT_DP}" | \
bcftools view --threads ${THREADS} -i "MEAN(FMT/DP) >= ${MIN_SITE_DP} && MEAN(FMT/DP) <= ${MAX_SITE_DP}" -Ou | \
bcftools +fill-tags -Ou -- -t AC,AN | \
bcftools view --threads ${THREADS} -i "AC > 0 && AC < AN && F_MISSING <= ${MAXMISS}" -Oz -o ${PREFIX}.snps.vcf.gz


bcftools index ${PREFIX}.snps.vcf.gz
echo "SNPs after filtering: $(bcftools index -n ${PREFIX}.snps.vcf.gz)"

# convert to phylip, heterozygotes as IUPAC codes, invariant columns removed
bcftools query -H -f '%REF\t%ALT[\t%GT]\n' ${PREFIX}.snps.vcf.gz | \
python3 ${SCRIPTDIR}/5a_vcf2phylip_snps.py ${PREFIX}.phy

# check the alignment
raxml-ng --check --msa ${PREFIX}.phy --model GTR+G+ASC_LEWIS --prefix ${PREFIX}.check

# compress to binary format for faster loading
raxml-ng --parse --msa ${PREFIX}.phy --model GTR+G+ASC_LEWIS --prefix ${PREFIX}

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
Rscript ${SCRIPTDIR}/5b_raxml_plot.r ${OUTDIR}/${PREFIX}.raxml.support ${SAMPLEINFO} ${OUTDIR} ${PREFIX} "${OUTGROUP}"
