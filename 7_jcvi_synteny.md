

```
micromamba create -n jcvi -c conda-forge -c bioconda jcvi minimap2 python=3.11
micromamba activate jcvi
micromamba install -c conda-forge ete4      # if conda-forge doesn't have it: pip install ete4
micromamba install -c conda-forge brotli      # if conda-forge doesn't have it: pip install ete4


PROJDIR=/nas/bcch_for_danny
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex

mkdir -p /nas/bcch_for_danny/analyses/synteny
bash ${SCRIPTDIR}/7_jcvi_synteny.sh

nohup bash ${SCRIPTDIR}/7_jcvi_synteny.sh > /nas/bcch_for_danny/analyses/synteny/synteny.log 2>&1 &
```