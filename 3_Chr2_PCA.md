# PCA of chromosome 2
A quick script to view clustering of autosomal diversity between species

Set up environment:
```
PROJDIR=/nas/bcch_for_danny/

SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex

REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna

CHR=Chr2
VCF=${PROJDIR}/data/vcfs/${CHR}.all_ind.diploid.vcf.gz
OUTDIR=${PROJDIR}/analyses/pca/${CHR}
SAMPLEINFO=${PROJDIR}/reference_lists/sample_info.txt
SAMPLES=${PROJDIR}/reference_lists/all_samples.txt
THREADS=8

mkdir -p ${OUTDIR}
cd ${OUTDIR}
```

Check that samples are assigned species etc.
```
# report samples in the vcf that are missing from the species lists, and vice versa
bcftools query -l ${VCF} | sort > vcf_samples.txt
sort ${SAMPLES} > list_samples.sorted.txt
echo "In vcf but not in species lists (will be dropped):"
comm -23 vcf_samples.txt list_samples.sorted.txt
echo "In species lists but not in vcf:"
comm -13 vcf_samples.txt list_samples.sorted.txt

# keep listed samples, biallelic snps, recompute AC/AN, drop sites with maf < 0.05
bcftools view --threads ${THREADS} -S ${SAMPLES} --force-samples ${VCF} -Ou | \
bcftools view --threads ${THREADS} -m2 -M2 -v snps -Ou | \
bcftools +fill-tags -Ou -- -t AF,AC,AN | \
bcftools view --threads ${THREADS} -q 0.05:minor -Oz -o ${CHR}.biallelic.maf05.vcf.gz

bcftools index ${CHR}.biallelic.maf05.vcf.gz

# convert to plink, assign unique snp ids, drop sites with > 20% missing data
plink --vcf ${CHR}.biallelic.maf05.vcf.gz \
--double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--snps-only 'just-acgt' \
--geno 0.2 \
--make-bed --out ${CHR}

# prune linked snps: 50kb window, step of 1 snp, r2 threshold of 0.5
plink --bfile ${CHR} --allow-extra-chr \
--indep-pairwise 50kb 1 0.5 \
--out ${CHR}_filtered

plink --bfile ${CHR} --allow-extra-chr \
--extract ${CHR}_filtered.prune.in \
--make-bed --out ${CHR}_pruned

echo "Variants before pruning: $(wc -l < ${CHR}.bim)"
echo "Variants after pruning: $(wc -l < ${CHR}_pruned.bim)"

# pca
plink --bfile ${CHR}_pruned --allow-extra-chr \
--pca 10 \
--out ${CHR}_pruned

# plot, colored by species
Rscript ${SCRIPTDIR}/3a_Chr2_PCA_plot.r ${OUTDIR}/${CHR}_pruned ${SAMPLEINFO} ${OUTDIR} ${CHR}
