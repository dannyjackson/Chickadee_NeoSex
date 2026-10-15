
micromamba create -n pretext -c conda-forge -c bioconda pretextmap samtools
micromamba activate pretext

samtools view -h /nas/Mountain.chickadee.pac.bio.reference/ready_reference/UNI3423.bam | PretextMap -o UNI3423.pretext --sortby length --mapq 10