#!/bin/bash
# B1_admixture.sh
# ADMIXTURE on the LD-pruned Chr2 plink files made by B1_PCA.sh, K = 2 to 6, with cross-validation
# usage: sbatch B1_admixture.sh   (or bash B1_admixture.sh)
# run B1_PCA.sh first

#SBATCH --job-name=admix_chr2
#SBATCH --output=slurm_output/admix_chr2.%j
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --time=24:00:00

# activate the micromamba env with admixture and R (ggplot2, dplyr, tidyr)
# install with: micromamba install -c conda-forge -c bioconda admixture r-base r-ggplot2 r-dplyr r-tidyr
ENV_NAME=your_env_name   # set to your env
eval "$(micromamba shell hook --shell bash)"
micromamba activate ${ENV_NAME}

if [ -z "${PROJDIR}" ]; then
    echo "Error: PROJDIR is not set" >&2
    exit 1
fi

CHR=Chr2
PCADIR=${PROJDIR}/analyses/pca/${CHR}
OUTDIR=${PROJDIR}/analyses/admixture/${CHR}
SAMPLEINFO=${PROJDIR}/reference_lists/sample_info.txt
SCRIPTDIR=$(dirname "$(readlink -f "$0")")   # if submitting with sbatch, set this to the directory holding B1_admixture_plot.r
THREADS=8
KMIN=2
KMAX=6

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# admixture only accepts integer chromosome codes, so copy the pruned files and recode the chromosome column
cp ${PCADIR}/${CHR}_pruned.bed ${CHR}_admix.bed
cp ${PCADIR}/${CHR}_pruned.fam ${CHR}_admix.fam
awk 'BEGIN{OFS="\t"} {$1=1; print}' ${PCADIR}/${CHR}_pruned.bim > ${CHR}_admix.bim

# run admixture for each K with 10-fold cross-validation
for K in $(seq ${KMIN} ${KMAX}); do
    admixture --cv=10 -j${THREADS} -s 12345 ${CHR}_admix.bed ${K} | tee ${CHR}_admix.K${K}.log
done

# collect cross-validation error, lowest is the best supported K
grep -h "CV error" ${CHR}_admix.K*.log | \
sed -E 's/CV error \(K=([0-9]+)\): ([0-9.]+)/\1\t\2/' | \
sort -n > ${CHR}_admix.cv_error.txt
echo "Cross-validation error by K:"
cat ${CHR}_admix.cv_error.txt

# plot, grouped by species
Rscript ${SCRIPTDIR}/B1_admixture_plot.r ${OUTDIR} ${CHR}_admix ${SAMPLEINFO} ${KMIN} ${KMAX} ${CHR}
