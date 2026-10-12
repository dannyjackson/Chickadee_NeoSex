# Managing chickadee environment
Install micromamba
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
Create environment for packages related to bcf/vcf management and pop gen analyses using genotype calls
```
export MAMBA_ROOT_PREFIX=$HOME/micromamba
echo 'export MAMBA_ROOT_PREFIX=$HOME/micromamba' >> ~/.bashrc

micromamba create -n bcf -c conda-forge -c bioconda bcftools samtools htslib
micromamba activate bcf
micromamba install -n bcf -c conda-forge -c bioconda plink
micromamba install -c conda-forge r-plotly r-htmlwidgets r-scatterplot3d pandoc
micromamba install -c conda-forge -c bioconda admixture r-tidyr
micromamba install -c conda-forge r-patchwork

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

micromamba install -n genespace -c conda-forge -c bioconda \
  "r-base>=4.5" bioconductor-biostrings bioconductor-rtracklayer \
  r-data.table r-igraph r-r.utils r-ggplot2 r-remotes r-rcpp


R
remotes::install_github("jtlovell/GENESPACE", dependencies = FALSE, upgrade = "never")


R.version.string    # confirm 4.5.x after the restore step
remotes::install_version("dbscan", version = "1.2.6",
                         repos = "https://cloud.r-project.org",
                         upgrade = "never")
packageVersion("dbscan")   # should say 1.2.6

library(GENESPACE)


micromamba create -n of255 -c conda-forge -c bioconda orthofinder=2.5.5 python=3.10 "numpy<2"
micromamba run -n of255 which orthofinder    # copy this path


micromamba env export -n genespace --explicit > ${PROJDIR}/reference_lists/genespace.lock.txt
micromamba env export -n of255 --explicit > ${PROJDIR}/reference_lists/of255.lock.txt

# environments can be recreated with lock files found under env_management github:
# micromamba create -n genespace -f genespace.lock.txt
# micromamba create -n of255 -f of255.lock.txt
```
PhaseWY
```
micromamba create -n Snakemake_PhaseWY -c conda-forge -c bioconda snakemake=8.27.1 conda
micromamba activate Snakemake_PhaseWY

cd /nas/bcch_for_danny/scripts/Snakemake_PhaseWY

cp -r slurm local
```
Repeat modeler
```
micromamba create -n repeats -c conda-forge -c bioconda repeatmodeler repeatmasker bedtools
micromamba activate repeats
RepeatModeler -version
```
# Managing project environment
Consistent variables across scripts
```
PROJDIR=/nas/bcch_for_danny/

SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex

REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
```