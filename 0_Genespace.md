# Genespace analysis

PROJDIR=/nas/bcch_for_danny/analyses/genespace/
SOURCEDIR=/nas/four_bird_synteny/genespace

mkdir -p input_files/gaga input_files/bcch input_files/moch input_files/zefi

cp $SOURCEDIR/bed/bcch.fa input_files/bcch/


for sp in bcch gaga moch zefi; do
  src=$sp
  [ "$sp" = "gaga" ] && src=chicken
  cp "$SOURCEDIR/bed/${src}.bed"    "input_files/${sp}/${sp}.bed"
  cp "$SOURCEDIR/peptide/${src}.fa" "input_files/${sp}/${sp}.fa"
done
ls input_files/*

micromamba activate genespace

```
cat > run_genespace.R << 'EOF'
library(GENESPACE)
gpar <- init_genespace(wd = "gs_run",
                       path2mcscanx = "~/MCScanX",
                       genomeIDs = c("gaga", "bcch", "moch", "zefi"),
                       nCores=12)
gpar <- run_genespace(gpar)
EOF

nohup Rscript run_genespace.R &
```