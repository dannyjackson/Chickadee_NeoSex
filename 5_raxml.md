#!/bin/bash
# B3_raxml.sh
# Maximum likelihood tree with raxml-ng from SNPs in one region of the original (unfiltered for PCA) Chr2 vcf
# region: CM059970.1:25440000-62170000
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
MAXMISS=0.2        # drop sites with more than this fraction of missing genotypes
BOOTSTRAPS=100
OUTGROUP=""        # optional, comma separated sample names to root the tree (e.g. BOCH samples); leave empty for midpoint rooting in the plot
PREFIX=${CHR}.${REGIONNAME}

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# region queries need an index on the input vcf
if [ ! -f ${VCF}.tbi ] && [ ! -f ${VCF}.csi ]; then
    echo "Error: no index for ${VCF}. Make one with: bcftools index -t ${VCF}" >&2
    exit 1
fi

# subset region and samples, keep biallelic snps that are variable in these samples, drop sites with too much missing data
bcftools view --threads ${THREADS} -r ${REGION} -S ${SAMPLES} --force-samples ${VCF} -Ou | \
bcftools view --threads ${THREADS} -m2 -M2 -v snps -Ou | \
bcftools +fill-tags -Ou -- -t AC,AN | \
bcftools view --threads ${THREADS} -i "AC > 0 && AC < AN && F_MISSING <= ${MAXMISS}" -Oz -o ${PREFIX}.snps.vcf.gz

bcftools index ${PREFIX}.snps.vcf.gz
echo "SNPs after filtering: $(bcftools index -n ${PREFIX}.snps.vcf.gz)"

# convert to phylip, heterozygotes as IUPAC codes, invariant columns removed
bcftools query -H -f '%REF\t%ALT[\t%GT]\n' ${PREFIX}.snps.vcf.gz | \
python3 ${SCRIPTDIR}/vcf2phylip_snps.py ${PREFIX}.phy

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
Rscript ${SCRIPTDIR}/B3_raxml_plot.r ${OUTDIR}/${PREFIX}.raxml.support ${SAMPLEINFO} ${OUTDIR} ${PREFIX} "${OUTGROUP}"
