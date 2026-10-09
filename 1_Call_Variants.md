# Generate the following VCFs:
1. All samples, aligned to raw ref genome, only call variants on Z.
    ***Use: compare depth of samples to confirm inferred sex***
2. All samples, aligned to ref genome with W removed, only call variants on Z using diploid calls
    ***Use: Phase variants on neo region of sex chromosome***
3. All females, aligned to raw ref genome, only call variants on W using haploid calls
    ***Use: Potentially extract W-linked variants in recent recombination suppressed region***


## Setup working environment
```
mkdir -p /nas/bcch_for_danny/data_from_will 
mv /nas/bcch_for_danny/* /nas/bcch_for_danny/data_from_will
mkdir -p /nas/bcch_for_danny/data /nas/bcch_for_danny/analyses /nas/bcch_for_danny/reference_lists /nas/bcch_for_danny/scripts
cd scripts/
git clone https://github.com/dannyjackson/Chickadee_NeoSex

micromamba activate bcf

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd ${PROJDIR}
REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna


```
# Call variants
## VCF 1: All samples, aligned to raw ref genome, only call variants on Z.
```
# make list of sorted bams
ls ${PROJDIR}/data_from_will/sex_linked_phylo/sorted_bam_files/*bam > ${PROJDIR}/reference_lists/sorted_bams.raw_genome.txt

mkdir -p ${PROJDIR}/data/vcfs/

# call variants on Z to confirm sex assignment
## Z chromosome
cat > scripts/call.Z.sh << 'EOF'
#!/bin/bash
micromamba activate bcf

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd ${PROJDIR}
REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna

bcftools mpileup -f ${REF} -r CM059970.1 -a AD,DP \
    -b ${PROJDIR}/reference_lists/sorted_bams.raw_genome.txt -Ou \
| bcftools call -mv --ploidy 2 -Oz -o ${PROJDIR}/data/vcfs/Z.all_ind.diploid.vcf.gz 

bcftools index ${PROJDIR}/data/vcfs/Z.all_ind.diploid.vcf.gz 

EOF

nohup bash scripts/call.Z.sh > scripts/call.Z.log 2>&1 &
```
## VCF 2: All samples, aligned to ref genome with W removed, only call variants on Z using diploid calls
```
# Remove W from ref genome
REF_noW=/nas/bcch_for_danny//data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.no_W.fna

samtools view ${REF} > ${REF_noW}

samtools faidx $REF

cut -f1 ${REF}.fai | grep -vxF CM059969.1 > $PROJDIR/reference_lists/scaffolds_noW.txt

samtools faidx $REF -r $PROJDIR/reference_lists/scaffolds_noW.txt > ${REF_noW}

samtools faidx ${REF_noW}


# make list of sorted bams
ls ${PROJDIR}/data_from_will/boch/sorted_bam_files/*bam > ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

ls ${PROJDIR}/data_from_will/moch/sorted_bam_files/*bam >> ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

ls ${PROJDIR}/data_from_will/cach/sorted_bam_files/*bam >> ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

# call variants on Z to confirm sex assignment
## Z chromosome
cat > scripts/call.Z_noW.sh << 'EOF'
#!/bin/bash
micromamba activate bcf

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd ${PROJDIR}
REF_noW=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.no_W.fna

bcftools mpileup -f ${REF_noW} -r CM059970.1 -a AD,DP \
    -b ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt -Ou \
| bcftools call -mv --ploidy 2 -Oz -o ${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz 

bcftools index ${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz 

EOF

nohup bash scripts/call.Z_noW.sh > scripts/call.Z_noW.log 2>&1 &