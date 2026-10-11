# Repeat modeler

micromamba activate repeats

PROJDIR=/nas/bcch_for_danny/

SCRIPTDIR=${PROJDIR}/scripts/

REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna

mkdir ${PROJDIR}/rep_mod && cd ${PROJDIR}/rep_mod

BuildDatabase -name bcch_ref_noW ${REF}

RepeatModeler -database chickadee -threads 32 -LTRStruct > rmod.log 2>&1
