# Infer the sex of each sample using depth of the Z

Setup environment
```

```

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
cd $SCRIPTDIR
git pull
nohup bash ${SCRIPTDIR}/2a_infer_sex.zdepth.sh &
```