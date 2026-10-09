# Managing chickadee environment
Install micromamba
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
Create environment for packages related to bcf/vcf management
```
micromamba create -n bcf -c conda-forge -c bioconda bcftools samtools htslib
micromamba activate bcf
```