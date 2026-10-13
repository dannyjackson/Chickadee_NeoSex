# Infer the sex of each sample using depth of the Z

Setup environment
```

```

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd $SCRIPTDIR
git pull

# bash ${SCRIPTDIR}/2a_infer_sex.zdepth.sh 
cd ${PROJDIR}

for sp in bcch cach moch boch; do
  mkdir -p ${PROJDIR}/analyses/sex_inference/${sp}
  nohup bash 2a_infer_sex.zdepth.sh \
    ${PROJDIR}/data/vcfs/Z_noW.${sp}.diploid.vcf.gz \
    ${PROJDIR}/data/vcfs/autosomes.all_ind.diploid.vcf.gz \
    > ${PROJDIR}/analyses/sex_inference/${sp}/infer_sex.log 2>&1 &
done

```