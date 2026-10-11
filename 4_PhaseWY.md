# PhaseWY 



```
micromamba activate Snakemake_PhaseWY

PROJDIR=/nas/bcch_for_danny/

SCRIPTDIR=${PROJDIR}/scripts/

REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.no_W.fna

VCF_BCCH=${PROJDIR}/data/vcfs/Z_noW.bcch.diploid.vcf.gz 
VCF_MOCH=${PROJDIR}/data/vcfs/Z_noW.moch.diploid.vcf.gz 

cd ${PROJDIR}
```
Make bed file for masking
```
python3 - ${REF} > ${PROJDIR}/reference_lists/repeats.bed <<'EOF'
import sys, re
name, seq = None, []
def emit(name, s):
    for m in re.finditer(r'[acgtn]+', s):
        print(f"{name}\t{m.start()}\t{m.end()}")
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            if name: emit(name, ''.join(seq))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name: emit(name, ''.join(seq))
EOF

sort -k1,1 -k2,2n ${PROJDIR}/reference_lists/repeats.bed > ${PROJDIR}/reference_lists/repeats.sorted.bed

/nas/bcch_for_danny/reference_lists/repeats.sorted.bed
```
Make sample info files
```
export to: /nas/bcch_for_danny/reference_lists/moch_samples.phaseWY_table.txt
In format:
sample_name depth   sex ploidy  bamfile
TT95871 34.8395 HETGAM  2   example/TT95871.bam
TT95873 26.9165 HOMGAM  2   example/TT95873.bam
```
snakemake -s ${SCRIPTDIR}/Snakemake_PhaseWY/workflow/PhaseWY --profile local --configfile ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml \
  --conda-prefix ~/snakemake_envs/phasewy -n   # dry run

snakemake -s${SCRIPTDIR}/Snakemake_PhaseWY/workflow/PhaseWY --profile local --configfile ${SCRIPTDIR}/Snakemake_PhaseWY/config/moch.yaml \
  --conda-prefix ~/snakemake_envs/phasewy