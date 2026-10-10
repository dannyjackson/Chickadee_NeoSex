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
library(GENESPACE)


gpar <- init_genespace(wd = "gs_run",
                       path2orthofinder = "/home/daja5529/micromamba/envs/of255/bin/orthofinder",
                       path2mcscanx = "/home/daja5529/micromamba/envs/genespace/bin",
                       genomeIDs = c("gaga", "bcch", "moch", "zefi"),
                       nCores = 1)

gpar <- run_genespace(gpar)

roi <- data.frame(
    genome = c("zefi", "zefi"),
    chr = c("Z", "1A"),
    color = c("#7f00e6", "#045b1b")
)

# Reload the run parameters if you're in a fresh R session
load("gs_run/results/gsParams.rda", verbose = TRUE)   # loads `gsParam`


pdf("gs_run/riparian/bcch_chrZ_syntenic.pdf", width = 10, height = 5)
ripZ <- plot_riparian(
  gsParam         = gsParam,
  refGenome       = "bcch",
  genomeIDs       = c("gaga", "bcch", "moch", "zefi"),
  highlightBed    = roi,
  useRegions      = TRUE,      # plot only regions syntenic to the ROI
  backgroundColor = NULL,      # drop all non-Z blocks
  useOrder        = FALSE      # x-axis in bp, not gene order
)
dev.off()


```