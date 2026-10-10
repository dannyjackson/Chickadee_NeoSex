# Scripts to make reference lists used across analyses

## General files
List of samples by species
```
cp ${PROJDIR}/data_from_will/moch/sample_list.txt  ${PROJDIR}/reference_lists/moch_samples.txt

cp ${PROJDIR}/data_from_will/new_data/sample_list.txt/sample_list.txt  ${PROJDIR}/reference_lists/bcch_samples.txt

cp ${PROJDIR}/data_from_will/cach/sample_list.txt  ${PROJDIR}/reference_lists/cach_samples.txt

```

## Files related to bcch ref genome
List of autosomes in raw bcch genome
```
awk '{print $1}' ${REF}.fai | grep 'CM' | grep -v '	CM059969.1' | grep -v 'CM059970.1' > ${PROJDIR}/reference_lists/bcch_ref.autosomes.txt
```

List of scaffolds in bcch ref genome excluding W
```
cut -f1 ${REF}.fai | grep -vxF CM059969.1 > $PROJDIR/reference_lists/scaffolds_noW.txt
```

## Bam files
List of sorted bams using alignment to raw bcch genome
```
ls ${PROJDIR}/data_from_will/sex_linked_phylo/sorted_bam_files/*bam > ${PROJDIR}/reference_lists/sorted_bams.raw_genome.txt
```

List of sorted bam files using alignment to bcch genome without W
```
ls ${PROJDIR}/data_from_will/boch/sorted_bam_files/*bam > ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

ls ${PROJDIR}/data_from_will/moch/sorted_bam_files/*bam >> ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

ls ${PROJDIR}/data_from_will/cach/sorted_bam_files/*bam >> ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

ls ${PROJDIR}/data_from_will/new_data/sorted_bam_files/*bam >> ${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt

```