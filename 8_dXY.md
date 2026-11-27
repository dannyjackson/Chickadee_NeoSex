# dXY across the Z: testing introgression against species-tree models

Asks whether BCCH and MOCH are most similar on the Z (and in the neo non-PAR) because:
- **Model 1:** introgression after the BCCH–CACH split. BCCH–MOCH Z divergence is *younger* than autosomal BCCH–CACH divergence (ratio < 1).
- **Model 2:** BCCH–MOCH is the species tree, and BCCH–CACH gene flow swamped the autosomes. BCCH–MOCH Z divergence is *older* than autosomal BCCH–CACH divergence (ratio > 1).

Both sides are scaled by divergence to BOCH, the outgroup, so the faster Z mutation rate (male mutation bias) cancels out.

The existing VCFs (`1a_Call_Variants.md`) are variant-only (`bcftools call -mv`), so 8a makes new all-sites calls. These include invariant sites, which pixy needs for the dXY denominator. Calls are on males only: ZZ males give true diploid Z genotypes, and the same males are used for the autosomes.

## Setup
```
micromamba activate bcf
micromamba install -n bcf -c conda-forge -c bioconda pixy pandas matplotlib

PROJDIR=/nas/bcch_for_danny
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
mkdir -p ${PROJDIR}/analyses/dxy
```

## 8a: all-sites calls and pixy (pi + dXY in 100 kb windows)
Run the Z first; it is much faster. Then run the autosomes.
```
cd ${PROJDIR}/analyses/dxy
DO_AUTOS=0 nohup bash ${SCRIPTDIR}/8a_dxy_allsites.sh > 8a.Z.log 2>&1 &
nohup bash ${SCRIPTDIR}/8a_dxy_allsites.sh > 8a.all.log 2>&1 &     # skips finished Z chunks
```
- Inputs:
  - male lists from `analyses/sex_inference/{sp}/sex_calls.tsv`
  - BAMs from `reference_lists/sorted_bams.noW_genome.txt`
  - the no-W reference
- Defaults:
  - 24 threads (`THREADS`)
  - 5 Mb calling chunks (`CHUNK`)
  - MAPQ and BASEQ ≥ 20
  - genotype DP ≥ 5 (`MIN_DP`)
  - site missingness ≤ 0.5 (`MAX_MISS`)
- Optional:
  - repeat mask: `NEGMASK=bed`
  - autosome subset: `AUTOS=list.txt`
- Restartable: finished chunks are skipped.

## 8b: compartments, outgroup scaling, block-bootstrap CIs, and the test
```
python3 ${SCRIPTDIR}/8b_dxy_summary.py pixy --out dxy_Z --mu 4.6e-9 --gen bcch=3.0,cach=2.3,moch=2.6
```
Z compartments use the PhaseWY / Genespace boundaries in BCCH coordinates:

| Compartment | Coordinates |
|---|---|
| neo-PAR | 0–14.78 and 25.44–62.17 Mb |
| MOCH-only stratum | 14.78–25.44 and 62.17–63.07 Mb |
| neo non-PAR | 63.07–65.20 Mb |
| ancestral non-PAR | 65.50–146.35 Mb; 65.50–70.94 Mb is reported separately because that stretch is not placed on the CACH Z |
| ancestral PAR | from 146.35 Mb to the end |

Windows that cross a boundary are dropped from the totals. "Autosomes (low π)" is the bottom quartile of autosomal windows by mean within-species π, used as a low-recombination proxy.

Outputs:

| File | Contents |
|---|---|
| `dxy_Z.windows.tsv` | per window: raw and outgroup-scaled dXY |
| `dxy_Z.summary.tsv` | per compartment × pair: dXY, scaled dXY, scaled Da, years, 95% CIs |
| `dxy_Z.test.tsv` | scaled Z BCCH–MOCH ÷ scaled autosomal BCCH–CACH, with the verdict |
| `dxy_Z.png` | (A) scaled dXY along the Z; (B) compartment totals |

Notes:
- **The test applies only to non-PAR compartments.** PAR compartments get a ratio for reference but no verdict.
- **The neo non-PAR is only ~2 Mb**, too small for 1 Mb blocks. Its CI is resampled by window, so it is too narrow, and the output flags it.
- **Smaller ancestral Ne on the Z (3/4)** lowers Z dXY under either model, slightly favouring a ratio < 1. Read the Da ratio alongside the dXY ratio, since Da is less affected.
- **Years scale with mu, which is still a placeholder.** The ratios do not depend on mu.
