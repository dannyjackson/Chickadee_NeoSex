# Subset variant call files
```
micromamba activate bcf

VCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz 

cd /nas/bcch_for_danny/data/vcfs

SAMPLES_BCCH=${PROJDIR}/reference_lists/bcch_samples.txt
SAMPLES_MOCH=${PROJDIR}/reference_lists/moch_samples.txt
SAMPLES_CACH=${PROJDIR}/reference_lists/cach_samples.txt
SAMPLES_BOCH=${PROJDIR}/reference_lists/boch_samples.txt

bcftools view ${VCF} -S ${SAMPLES_BCCH} | bgzip > ${PROJDIR}/data/vcfs/Z_noW.bcch.diploid.vcf.gz 

bcftools view ${VCF} -S ${SAMPLES_MOCH} | bgzip > ${PROJDIR}/data/vcfs/Z_noW.moch.diploid.vcf.gz 

bcftools view ${VCF} -S ${SAMPLES_CACH} | bgzip > ${PROJDIR}/data/vcfs/Z_noW.cach.diploid.vcf.gz 

bcftools view ${VCF} -S ${SAMPLES_BOCH} | bgzip > ${PROJDIR}/data/vcfs/Z_noW.boch.diploid.vcf.gz 

```