# MSMC2 + RCCR from haplotype-resolved assemblies (BCCH, CACH, MOCH)
Adapted from the Sulidae MSMC2 → RCCR workflow (SulidaeGenomics/Analyses/3_Demographic_History/3c_MSMC_Divergence/MSMC.md and RCCR_BFPE_MANA_Temperature.*.R), with these changes:

| Step | Sulidae (short reads, several individuals per species) | Chickadees (hap1 + hap2 assemblies, one individual per species) |
|---|---|---|
| Genotypes | bcftools + bamCaller.py per individual | dipcall: align hap1 and hap2 to the shared BCCH reference |
| Callability mask | individual depth mask + SNPable mask | dipcall confident-region BED (+ optional SNPable / repeat masks), split per chromosome |
| Phasing | WhatsHap | Already phased: GT = hap1\|hap2 |
| Sex chromosomes | drop W/Z (NC_087547.1 / NC_087548.1) | drop W (CM059969.1) and Z incl. neo-Z (CM059970.1) |
| Input files | one multihetsep per species per chromosome | **one joint multihetsep per chromosome containing all species** |
| Within-species MSMC2 | `-I 0,1,...,2n-1` on that species' files | `-I 0,1` / `2,3` / `4,5` on the joint files |
| Cross MSMC2 | `-I 0,0,...,1,1,...` on BFBO + PEBO files | `-s -I 0-2,0-3,1-2,1-3` (etc.) on the joint files |
| Time pattern / iterations | `-p 1*2+15*1+1*2 -i 100` | same |
| Combine / RCCR | combineCrossCoal.py → RCCR, 3 generation times, LGM/PGM | same (6e mirrors the Sulidae RCCR figure) |
| New (optional) | – | bootstraps (6d), MSMC-IM migration (RUN_IM=1) |

**Why the cross step differs:** msmc2 reads a plain list for `-I` (e.g. `0,0,0,1,1`) as "every pair of the listed positions". In the Sulidae runs that produced pairs 0–0, 0–1 and 1–1 only: one haplotype compared with itself, or the two haplotypes of the first individual in each file. Because the BFBO and PEBO files were separate, no pair crossed species. Cross-coalescence needs both species' haplotypes in the same multihetsep file, plus explicit pairs that cross species (`-I a-c,a-d,b-c,b-d`).

All species are aligned to the same reference (bPoeAtr1.hap1, GCA_030490865.1), so the joint input files share coordinates. Haplotype indices: **bcch = 0,1; cach = 2,3; moch = 4,5**.

## Setup environment
```
micromamba create -n msmc -c conda-forge -c bioconda msmc2 bcftools htslib samtools python=3.11 numpy scipy matplotlib
micromamba activate msmc

mkdir -p ~/programs && cd ~/programs
git clone https://github.com/stschiff/msmc-tools
git clone https://github.com/wangke16/MSMC-IM
wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2
tar -jxf dipcall-0.3_x64-linux.tar.bz2      # creates ~/programs/dipcall.kit (bundles minimap2, k8, htsbox, bedtk)
```

## Consistent variables
```
PROJDIR=/nas/bcch_for_danny
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
BCCH_HAP1=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
BCCH_HAP2=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490855.1/GCA_030490855.1_bPoeAtr1.hap2_genomic.fna
CACH_HAP1=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap1.p_ctg.fa.gz
CACH_HAP2=/nas/bcch_for_danny/cach_reads/bPoeCar1.hap2.p_ctg.fa.gz
MOCH_HAP1=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/Poecile_gambeli_hap1_NCBI.fasta
MOCH_HAP2=/nas/Mountain.chickadee.pac.bio.reference/ready_reference/ncbi_hap2_submission/aspera_upload/Poecile_gambeli_hap2.fasta
```
(The haplotype paths are also hard-coded in `6a_dipcall_haps.sh`.)

## Download BCCH hap2
```
FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

cd ${PROJDIR}/data_from_will
datasets download genome accession GCA_030490855.1 \
  --include "${FILES_TO_DOWNLOAD}" \
  --filename GCA_030490855.1.zip

unzip GCA_030490855.1.zip
```
## 1. Call phased variants from each pair of haplotypes (dipcall)
```
micromamba activate msmc
for sp in bcch cach moch; do
  mkdir -p ${PROJDIR}/analyses/msmc_rccr/dipcall/${sp}
  nohup bash ${SCRIPTDIR}/6a_dipcall_haps.sh ${sp} \
    > ${PROJDIR}/analyses/msmc_rccr/dipcall/${sp}/dipcall.log 2>&1 &
done
```
Check before moving on:
```
for sp in bcch cach moch; do
  echo ${sp}; tail -2 ${PROJDIR}/analyses/msmc_rccr/dipcall/${sp}/dipcall.log
done
```
- **Callable fraction** should be similar across species. If MOCH is much lower than BCCH/CACH, the divergent alignments are failing under `asm5`: delete `dipcall/moch/` and rerun with `ASM_PRESET=asm20 bash ${SCRIPTDIR}/6a_dipcall_haps.sh moch`.
- **Heterozygosity** should be roughly comparable to read-based estimates for these species. Values near 0 mean hap1 = hap2 (wrong file), and very high values usually mean collapsed or duplicated sequence.

## 2. Build joint MSMC2 input (all three species, per autosome) using ref genome's soft mask
```
PROJDIR=/nas/bcch_for_danny
REF=${PROJDIR}/data_from_will/ncbi_dataset/data/GCA_030490865.1/GCA_030490865.1_bPoeAtr1.hap1_genomic.fna
python3 ${SCRIPTDIR}/softmask2bed.py ${REF} > ${PROJDIR}/analyses/msmc_rccr/bPoeAtr1.hap1.softmask.bed

NEGMASK=${PROJDIR}/analyses/msmc_rccr/bPoeAtr1.hap1.softmask.bed \
  nohup bash ${SCRIPTDIR}/6b_msmc_input.sh > ${PROJDIR}/analyses/msmc_rccr/msmc_input.log 2>&1 &
  
```

## 3. Run MSMC2 within and across species, combine, MSMC-IM
Defaults match the Sulidae runs (`-p 1*2+15*1+1*2 -i 100`). Add `RUN_IM=1 MU=<rate>` to also fit MSMC-IM; `MU` must equal `mu` in `6e`.
```
MSMC2=msmc2_Linux THREADS=24 nohup bash ${SCRIPTDIR}/6c_run_msmc_rccr.sh \
  > ${PROJDIR}/analyses/msmc_rccr/run_msmc.log 2>&1 &
```
Outputs in `analyses/msmc_rccr/output/`:
- `within.{bcch,cach,moch}.final.txt`: Ne through time
- `cross.{pair}.final.txt`: cross-coalescence rate
- `combined.{pair}.final.txt`: lambda_00 / lambda_01 / lambda_11 → RCCR
- `msmcim/{pair}.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt`: m(t), M(t) (only with RUN_IM=1; time in generations)

## 4. Bootstraps (optional; not in the Sulidae workflow)
Each replicate reruns all six MSMC2 runs, so start with a few:
```
BOOT_IDS="1 2 3 4 5" THREADS=24 nohup bash ${SCRIPTDIR}/6d_msmc_bootstraps.sh \
  > ${PROJDIR}/analyses/msmc_rccr/bootstraps.log 2>&1 &
```

## 5. Plot locally
```
scp -r daja5529@chickadee.colorado.edu:/nas/bcch_for_danny/analyses/msmc_rccr/output ~/Desktop/chickadee_msmc_rccr/
scp -r daja5529@chickadee.colorado.edu:/nas/bcch_for_danny/analyses/msmc_rccr/bootstrap/outputs ~/Desktop/chickadee_msmc_rccr/bootstrap/
```
Then run `6e_plot_msmc_rccr.r` in RStudio. **Set `mu` and the three generation times `gens` at the top first; the defaults are placeholders.** To get the temperature panel, put a copy of the Sulidae `temperatures.csv` in the output folder.

## Caveats
- With one individual (2 haplotypes) per species, MSMC2 has little resolution in the recent past. Interpret only the RCCR decline itself and anything older.
- Every species is mapped to the BCCH reference. Regions that fail to align in a divergent species are masked for all species (masks are intersected), which keeps within and cross rates comparable but reduces total sequence.
- CACH haplotypes are hifiasm contigs (`p_ctg`), and MOCH hap1/hap2 were scaffolded separately with HiRise. Phase can therefore switch between contigs. This has no effect on within-species runs and only a small effect on cross-coalescence.
