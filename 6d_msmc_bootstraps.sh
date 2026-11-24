#!/bin/bash
# -----------------------------------------------------------------------------
# Block bootstraps for the within, cross and combined MSMC2 runs.
# Resamples 3 Mb chunks of the joint input files (6b) into NBOOT pseudo-genomes,
# then reruns 6c_run_msmc_rccr.sh on each one (MSMC-IM off by default).
#
# Pseudo-genome size = NR_CHR x CHUNKS_PER_CHR x CHUNK bp. The defaults
# (18 x 20 x 3 Mb = 1.08 Gb) roughly match the callable chickadee autosomes.
#
# Usage:
#   bash 6d_msmc_bootstraps.sh            # generate + run all bootstraps
#   BOOT_IDS="1 2 3" bash 6d_msmc_bootstraps.sh   # run only some replicates
#
# Optional env overrides:
#   SPECIES="bcch cach moch"  NBOOT=20  CHUNK=3000000  CHUNKS_PER_CHR=20
#   NR_CHR=18  SEED=1  BOOT_IDS="1 ... NBOOT"  BOOT_IM=0
#   plus anything 6c_run_msmc_rccr.sh accepts (THREADS, P_PAR, ...)
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny
SCRIPTDIR=${SCRIPTDIR:-${PROJDIR}/scripts/Chickadee_NeoSex}
MSMCTOOLS=${MSMCTOOLS:-${HOME}/programs/msmc-tools}
SPECIES=${SPECIES:-"bcch cach moch"}
NBOOT=${NBOOT:-20}
CHUNK=${CHUNK:-3000000}
CHUNKS_PER_CHR=${CHUNKS_PER_CHR:-20}
NR_CHR=${NR_CHR:-18}
SEED=${SEED:-1}
BOOT_IDS=${BOOT_IDS:-$(seq 1 ${NBOOT})}
BOOT_IM=${BOOT_IM:-0}
OUTBASE=${PROJDIR}/analyses/msmc_rccr
TAG=$(echo ${SPECIES} | tr ' ' '_')
BSDIR=${OUTBASE}/bootstrap

mkdir -p ${BSDIR}

# --- 1. Generate bootstrapped inputs (once) ---
# out_dir_prefix and files are POSITIONAL arguments of multihetsep_bootstrap.py
if [[ ! -d ${BSDIR}/${TAG}.bootstrap_${NBOOT} ]]; then
  echo "[$(date)] Generating ${NBOOT} bootstraps"
  python3 ${MSMCTOOLS}/multihetsep_bootstrap.py \
    -n ${NBOOT} -s ${CHUNK} --chunks_per_chromosome ${CHUNKS_PER_CHR} \
    --nr_chromosomes ${NR_CHR} --seed ${SEED} \
    ${BSDIR}/${TAG}.bootstrap $(cat ${OUTBASE}/input/${TAG}.input_files.txt)
else
  echo "[$(date)] Bootstrap inputs exist in ${BSDIR}; reusing"
fi

# --- 2. Run within / cross / combined on each replicate ---
for b in ${BOOT_IDS}; do
  echo "[$(date)] Bootstrap ${b}"
  SPECIES="${SPECIES}" RUN_IM=${BOOT_IM} \
  INPUT_DIR=${BSDIR}/${TAG}.bootstrap_${b} \
  INPUT_GLOB="bootstrap_multihetsep.chr*.txt" \
  OUTPUT_DIR=${BSDIR}/outputs/boot_${b} \
    bash ${SCRIPTDIR}/6c_run_msmc_rccr.sh
done

echo "[$(date)] Finished bootstraps ${BOOT_IDS}"
