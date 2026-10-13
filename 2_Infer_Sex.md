# Infer the sex of each sample using depth of the Z

Setup environment
```
PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd $SCRIPTDIR
git pull
cd ${PROJDIR}
micromamba activate bcf
```
Infer sex from VCF
```
for sp in bcch cach moch boch; do
  mkdir -p ${PROJDIR}/analyses/sex_inference/${sp}
  bash ${SCRIPTDIR}/2a_infer_sex.zdepth.vcf.sh \
    ${PROJDIR}/data/vcfs/Z_noW.${sp}.diploid.vcf.gz \
    ${PROJDIR}/data/vcfs/autosomes.all_ind.diploid.vcf.gz \
    > ${PROJDIR}/analyses/sex_inference/${sp}/infer_sex.log 2>&1 &
done
```
Infer sex from bam
```
for sp in bcch cach moch boch; do
  mkdir -p ${PROJDIR}/analyses/sex_inference/${sp}
  nohup bash ${SCRIPTDIR}/2a_infer_sex.zdepth.bam.sh ${sp} \
    > ${PROJDIR}/analyses/sex_inference/${sp}/infer_sex.log 2>&1 &
done
```