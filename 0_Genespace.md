# Genespace analysis
Setup environment
```
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
Run orthofinder independently
```
micromamba run -n of255 python --version    # should say 3.10.x
micromamba run -n of255 orthofinder -f gs_run/tmp -t 12 -a 1 -X -o gs_run/orthofinder
```
Run genespace
```
cat > run_genespace.R << 'EOF'
library(GENESPACE)
gpar <- init_genespace(wd = "gs_run",
                       path2orthofinder = "/home/daja5529/micromamba/envs/of255/bin/orthofinder",
                       path2mcscanx = "/home/daja5529/MCScanX",
                       genomeIDs = c("gaga", "bcch", "moch", "zefi"),
                       nCores=12)
gpar <- run_genespace(gpar)
EOF

nohup Rscript run_genespace.R &
```