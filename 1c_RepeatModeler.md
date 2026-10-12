# Repeat modeler

micromamba activate repeats

PROJDIR=/nas/bcch_for_danny/

SCRIPTDIR=${PROJDIR}/scripts/

REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna

mkdir ${PROJDIR}/rep_mod 


mkdir -p /tmp/rmod_bcch && cd /tmp/rmod_bcch
BuildDatabase -name bcch_ref "$REF"

RepeatModeler -database bcch_ref -threads 32 -LTRStruct > rmod.log 2>&1


&& cd ${PROJDIR}/rep_mod