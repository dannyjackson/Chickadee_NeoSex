# Scripts to make reference lists used across analyses

## General files
List of samples by species
```
cp ${PROJDIR}/data_from_will/moch/sample_list.txt  ${PROJDIR}/reference_lists/moch_samples.txt

cp ${PROJDIR}/data_from_will/new_data/sample_list.txt  ${PROJDIR}/reference_lists/bcch_samples.txt

cp ${PROJDIR}/data_from_will/cach/sample_list.txt  ${PROJDIR}/reference_lists/cach_samples.txt

ls ../data_from_will/boch/bam_files/ | awk -F "." '{print $1}' >  ${PROJDIR}/reference_lists/boch_samples.txt



# combine into one tab-delimited file: Sample  Species
# takes the first column of each list, strips windows line endings and blank lines
SAMPLEINFO=${PROJDIR}/reference_lists/sample_info.txt
echo -e "Sample\tSpecies" > ${SAMPLEINFO}

for sp in moch bcch cach boch; do
    SPUP=$(echo ${sp} | tr '[:lower:]' '[:upper:]')
    tr -d '\r' < ${PROJDIR}/reference_lists/${sp}_samples.txt | \
    awk -v sp=${SPUP} 'NF > 0 {print $1"\t"sp}' >> ${SAMPLEINFO}
done

# report counts per species
echo "Samples per species:"
tail -n +2 ${SAMPLEINFO} | cut -f2 | sort | uniq -c

# flag any sample assigned to more than one species
DUPS=$(tail -n +2 ${SAMPLEINFO} | cut -f1 | sort | uniq -d)
if [ -n "${DUPS}" ]; then
    echo "Warning: these samples appear in more than one species list:" >&2
    echo "${DUPS}" >&2
fi

# list of all samples, used to subset the vcf
tail -n +2 ${SAMPLEINFO} | cut -f1 > ${PROJDIR}/reference_lists/all_samples.txt
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