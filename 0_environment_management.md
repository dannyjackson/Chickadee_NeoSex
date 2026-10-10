# Managing chickadee environment
Install micromamba
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
Create environment for packages related to bcf/vcf management
```
export MAMBA_ROOT_PREFIX=$HOME/micromamba
echo 'export MAMBA_ROOT_PREFIX=$HOME/micromamba' >> ~/.bashrc

micromamba create -n bcf -c conda-forge -c bioconda bcftools samtools htslib
micromamba activate bcf
```
Create environment for genespace
```

micromamba create -n genespace -c conda-forge -c bioconda r-base=4.3 orthofinder=2.5.5 r-devtools r-biocmanager
micromamba activate genespace

micromamba install -n genespace -c conda-forge -c bioconda \
  r-remotes r-data.table r-igraph r-dbscan r-r.utils r-ggplot2 \
  bioconductor-biostrings bioconductor-rtracklayer

micromamba install -n genespace -c conda-forge -c bioconda orthofinder=2.5.5

micromamba install -n genespace -c conda-forge -c bioconda mcscanx

micromamba install -n genespace -c conda-forge "r-dbscan<1.2"

R
remotes::install_github("jtlovell/GENESPACE", dependencies = FALSE, upgrade = "never")
library(GENESPACE)



micromamba create -n of255 -c conda-forge -c bioconda orthofinder=2.5.5 python=3.10 "numpy<2"
micromamba run -n of255 which orthofinder    # copy this path
```