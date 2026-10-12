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
Compute average depth from autosomal vcf
```
ZVCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz
OUT=${PROJDIR}/data/vcfs/depth_check        # new folder; nothing existing is touched
mkdir -p "$OUT"

# Mean depth per sample
bcftools view -h "$ZVCF" | grep -m1 'ID=DP,.*FORMAT\|FORMAT=<ID=DP' || echo "No FORMAT/DP in header!"

vcftools --gzvcf "$ZVCF" --depth --out "$OUT/zchr"
# -> $OUT/zchr.idepth : INDV  N_SITES  MEAN_DEPTH
sort -k3,3n "$OUT/zchr.idepth" | column -t
```
Make sample info files
```
export to: /nas/bcch_for_danny/reference_lists/moch_samples.phaseWY_table.txt
In format:
sample_name depth   sex ploidy  bamfile
TT95871 34.8395 HETGAM  2   example/TT95871.bam
TT95873 26.9165 HOMGAM  2   example/TT95873.bam

################################################
# file for moch
################################################

SEX=/nas/bcch_for_danny/analyses/sex_inference/sex_calls.tsv
BAMS=/nas/bcch_for_danny/reference_lists/sorted_bams.noW_genome.txt
DEPTH=/nas/bcch_for_danny/data/vcfs/depth_check/zchr.idepth
OUT=/nas/bcch_for_danny/reference_lists/moch_samples.phaseWY_table.txt
SPECIES=moch

# 1. Build the table
{
  printf "sample_name\tdepth\tsex\tploidy\tbamfile\n"
  awk -v OFS='\t' '
    FILENAME==ARGV[1] { if (FNR>1) sex[$1]=$NF; next }   # sex = last column
    FILENAME==ARGV[2] { if (FNR>1) dp[$1]=$3;  next }    # MEAN_DEPTH
    {
      bam=$0; n=split(bam,p,"/"); s=p[n]; sub(/_sorted.*/,"",s)
      if (!(s in sex)) { print "SKIP (no sex call): " s > "/dev/stderr"; next }
      if (!(s in dp))  { print "SKIP (no depth): "    s > "/dev/stderr"; next }
      if (sex[s]=="F") g="HETGAM"; else if (sex[s]=="M") g="HOMGAM"
      else { print "SKIP (sex=" sex[s] "): " s > "/dev/stderr"; next }
      print s, dp[s], g, 2, bam
    }' "$SEX" "$DEPTH" <(grep "/${SPECIES}/" "$BAMS")
} > "$OUT"

# 2. Summary
column -t "$OUT"
echo "Females: $(grep -c HETGAM "$OUT")   Males: $(grep -c HOMGAM "$OUT")"


################################################
# File for bcch
################################################
OUT=/nas/bcch_for_danny/reference_lists/bcch_samples.phaseWY_table.txt
SPECIES=new_data

# 1. Build the table
{
  printf "sample_name\tdepth\tsex\tploidy\tbamfile\n"
  awk -v OFS='\t' '
    FILENAME==ARGV[1] { if (FNR>1) sex[$1]=$NF; next }   # sex = last column
    FILENAME==ARGV[2] { if (FNR>1) dp[$1]=$3;  next }    # MEAN_DEPTH
    {
      bam=$0; n=split(bam,p,"/"); s=p[n]; sub(/_sorted.*/,"",s)
      if (!(s in sex)) { print "SKIP (no sex call): " s > "/dev/stderr"; next }
      if (!(s in dp))  { print "SKIP (no depth): "    s > "/dev/stderr"; next }
      if (sex[s]=="F") g="HETGAM"; else if (sex[s]=="M") g="HOMGAM"
      else { print "SKIP (sex=" sex[s] "): " s > "/dev/stderr"; next }
      print s, dp[s], g, 2, bam
    }' "$SEX" "$DEPTH" <(grep "/${SPECIES}/" "$BAMS")
} > "$OUT"

# 2. Summary
column -t "$OUT"
echo "Females: $(grep -c HETGAM "$OUT")   Males: $(grep -c HOMGAM "$OUT")"
```
Make list of contigs in vcf
```
bcftools index -t /nas/bcch_for_danny/data/vcfs/Z_noW.moch.diploid.vcf.gz
bcftools index -s /nas/bcch_for_danny/data/vcfs/Z_noW.moch.diploid.vcf.gz | cut -f1 \
  > /nas/bcch_for_danny/reference_lists/moch_contigs.phaseWY.txt
cat /nas/bcch_for_danny/reference_lists/moch_contigs.phaseWY.txt   # should be the Z only

```
Filter vcf to sites that are only invariant in each species, and filter for mapping quality
```
####################################
# moch
####################################
IN=/nas/bcch_for_danny/data/vcfs/Z_noW.moch.diploid.vcf.gz
OUT=/nas/bcch_for_danny/data/vcfs/Z_noW.moch.diploid.biSNP.vcf.gz

bcftools view -a -Ou "$IN" \
 | bcftools view -m2 -M2 -v snps -i 'AC>0 && AC<AN && MQ>=30 && QUAL>=30' -Oz -o "$OUT"
bcftools index -t "$OUT"

bcftools view -H -i 'AC==0 || AC==AN' "$OUT" | head -3   # should print nothing

###################################
# bcch
###################################
IN=/nas/bcch_for_danny/data/vcfs/Z_noW.bcch.diploid.vcf.gz
OUT=/nas/bcch_for_danny/data/vcfs/Z_noW.bcch.diploid.biSNP.vcf.gz

bcftools view -a -Ou "$IN" \
 | bcftools view -m2 -M2 -v snps -i 'AC>0 && AC<AN && MQ>=30 && QUAL>=30' -Oz -o "$OUT"
bcftools index -t "$OUT"

bcftools view -H -i 'AC==0 || AC==AN' "$OUT" | head -3   # should print nothing
```
# dry run
snakemake -s ${SCRIPTDIR}/Snakemake_PhaseWY/workflow/PhaseWY --profile ${SCRIPTDIR}/Snakemake_PhaseWY/local --configfile ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml \
  --conda-prefix ~/snakemake_envs/phasewy -n   # dry run

# run with dist_model: "Hamming" and again with "Inverse_MAF"
snakemake -s${SCRIPTDIR}/Snakemake_PhaseWY/workflow/PhaseWY --profile ${SCRIPTDIR}/Snakemake_PhaseWY/local --configfile ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml \
  --conda-prefix ~/snakemake_envs/phasewy --rerun-incomplete

# make additional yamls
## make yaml for bcch
grep 'moch' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml
sed 's/moch/bcch/g' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml > ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.bcch.yaml

# make dist_model Inverse_MAF yaml for moch and bcch
grep 'Hamming' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml
grep 'Hamming' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.bcch.yaml

sed -i 's/Hamming/Inverse_MAF/g' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml
sed -i 's/Hamming/Inverse_MAF/g' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.bcch.yaml

# set back

sed -i 's/Inverse_MAF/Hamming/g' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.moch.yaml
sed -i 's/Inverse_MAF/Hamming/g' ${SCRIPTDIR}/Snakemake_PhaseWY/config/config.bcch.yaml


# Visualize snps on sex limted chr
cd ${SCRIPTDIR}/Chickadee_NeoSex
git pull

cd /nas/bcch_for_danny/scripts/Snakemake_PhaseWY/results/moch_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming
# 1. Pull position, allele count and genetic position (fast; skips genotypes)
micromamba activate bcf

bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/CM\n' vcfs/sexlimited.vcf.gz > sex_limited_sites.tsv

Rscript ${SCRIPTDIR}/Chickadee_NeoSex/4a_plot_z_snp_density.R sex_limited_sites.tsv CM059970.1 100000