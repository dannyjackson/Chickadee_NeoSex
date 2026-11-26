# PSMC

Define environment variables
```
PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
WORK
GCA_030490855.1

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

cd ${PROJDIR}/data_from_will
datasets download genome accession GCA_030490855.1 \
  --include "${FILES_TO_DOWNLOAD}" \
  --filename GCA_030490855.1.zip

unzip GCA_030490855.1.zip

BCCH_HAP1 = "${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna"
BCCH_HAP2 = "${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490855.1/GCA_030490855.1_bPoeAtr1.hap2_genomic.fna"
MOCH_HAP1="/nas/Mountain.chickadee.pac.bio.reference/ready_reference/Poecile_gambeli_hap1_NCBI.fasta"
MOCH_HAP2="/nas/Mountain.chickadee.pac.bio.reference/ready_reference/ncbi_hap2_submission/aspera_upload/Poecile_gambeli_hap2.fasta"
CACH_HAP1 = "/nas/bcch_for_danny/cach_reads/bPoeCar1.hap1.p_ctg.fa.gz"
CACH_HAP2 = "/nas/bcch_for_danny/cach_reads/bPoeCar1.hap2.p_ctg.fa.gz"
```

/nas/Mountain.chickadee.pac.bio.reference/ready_reference/UNI3423.bam