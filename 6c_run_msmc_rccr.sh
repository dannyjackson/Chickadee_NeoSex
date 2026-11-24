#!/bin/bash
# -----------------------------------------------------------------------------
# Run MSMC2 within each species and across each species pair, then combine
# into the files used for RCCR (and optionally MSMC-IM).
#
# Each species is one individual = 2 haplotypes:
#   within <sp>     : msmc2 -I a,b                       (Ne of that species)
#   cross  <sp1_sp2>: msmc2 -s -I a-c,a-d,b-c,b-d       (cross-coalescence only)
#   combined        : combineCrossCoal.py cross within1 within2
#                     -> columns lambda_00, lambda_01, lambda_11
#   RCCR = 2 * lambda_01 / (lambda_00 + lambda_11)       (computed in 6e plot)
#
# All runs use the same JOINT input files from 6b_msmc_input.sh (every species'
# haplotypes in one file per chromosome), so within and cross rates come from
# exactly the same callable sites. Cross-coalescence needs this: msmc2 can only
# pair haplotypes that sit in the same input file.
#
# Usage:
#   bash 6c_run_msmc_rccr.sh
#
# Optional env overrides:
#   SPECIES="bcch cach moch"  (must match the order used in 6b)
#   MSMC2=msmc2  MSMCTOOLS=~/programs/msmc-tools  MSMCIM=~/programs/MSMC-IM
#   THREADS=24  NUM_OPT=100  P_PAR="1*2+15*1+1*2"   (as in the Sulidae runs)
#   RUN_IM=0  MU=4.6e-9       (RUN_IM=1 adds MSMC-IM; MU per site per generation)
#   INPUT_DIR=...  OUTPUT_DIR=...   (used by 6d to rerun on bootstraps)
# Runs whose .final.txt already exists are skipped; delete to rerun.
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny
MSMC2=${MSMC2:-msmc2}
MSMCTOOLS=${MSMCTOOLS:-${HOME}/programs/msmc-tools}
MSMCIM=${MSMCIM:-${HOME}/programs/MSMC-IM}
SPECIES=${SPECIES:-"bcch cach moch"}
THREADS=${THREADS:-24}
NUM_OPT=${NUM_OPT:-100}
P_PAR=${P_PAR:-"1*2+15*1+1*2"}
RUN_IM=${RUN_IM:-0}
MU=${MU:-4.6e-9}
OUTBASE=${PROJDIR}/analyses/msmc_rccr
TAG=$(echo ${SPECIES} | tr ' ' '_')
INPUT_DIR=${INPUT_DIR:-${OUTBASE}/input}
OUTPUT_DIR=${OUTPUT_DIR:-${OUTBASE}/output}
INPUT_GLOB=${INPUT_GLOB:-msmc_input.${TAG}.*.txt}

mkdir -p ${OUTPUT_DIR}/log_and_loop

INPUTS=$(ls ${INPUT_DIR}/${INPUT_GLOB} 2>/dev/null || true)
[[ -n ${INPUTS} ]] || { echo "ERROR: no inputs matching ${INPUT_DIR}/${INPUT_GLOB}"; exit 1; }

# haplotype indices: species k -> 2k, 2k+1
declare -A H1 H2
SP_ARR=(${SPECIES})
for k in "${!SP_ARR[@]}"; do
  H1[${SP_ARR[$k]}]=$((2*k))
  H2[${SP_ARR[$k]}]=$((2*k+1))
done

date
echo "Script: 6c_run_msmc_rccr.sh"
echo "Species / haplotypes: $(for s in ${SPECIES}; do printf '%s=%s,%s ' $s ${H1[$s]} ${H2[$s]}; done)"
echo "Inputs: $(echo ${INPUTS} | wc -w) files in ${INPUT_DIR}"
echo "Time pattern: ${P_PAR}; iterations: ${NUM_OPT}; threads: ${THREADS}"
echo "Output dir: ${OUTPUT_DIR}"

run_msmc () {   # $1 = output prefix, $2 = index string, $3 = extra flags
  local OUT=${OUTPUT_DIR}/$1
  if [[ -s ${OUT}.final.txt ]]; then
    echo "[$(date)] $1 exists, skipping"; return
  fi
  echo "[$(date)] Running $1 (-I $2 $3)"
  ${MSMC2} -t ${THREADS} -p "${P_PAR}" -i ${NUM_OPT} $3 -I $2 -o ${OUT} ${INPUTS}
  mv ${OUT}.loop.txt ${OUT}.log ${OUTPUT_DIR}/log_and_loop/ 2>/dev/null || true
}

# --- 1. Within-species coalescence (Ne) ---
for sp in ${SPECIES}; do
  run_msmc within.${sp} "${H1[$sp]},${H2[$sp]}" ""
done

# --- 2. Cross-species coalescence, then combine ---
for ((a=0; a<${#SP_ARR[@]}; a++)); do
  for ((b=a+1; b<${#SP_ARR[@]}; b++)); do
    s1=${SP_ARR[$a]}; s2=${SP_ARR[$b]}; PAIR=${s1}_${s2}
    IDX="${H1[$s1]}-${H1[$s2]},${H1[$s1]}-${H2[$s2]},${H2[$s1]}-${H1[$s2]},${H2[$s1]}-${H2[$s2]}"
    run_msmc cross.${PAIR} "${IDX}" "-s"

    COMB=${OUTPUT_DIR}/combined.${PAIR}.final.txt
    if [[ ! -s ${COMB} ]]; then
      python3 ${MSMCTOOLS}/combineCrossCoal.py \
        ${OUTPUT_DIR}/cross.${PAIR}.final.txt \
        ${OUTPUT_DIR}/within.${s1}.final.txt \
        ${OUTPUT_DIR}/within.${s2}.final.txt > ${COMB}
      echo "[$(date)] Combined -> ${COMB}"
    fi

    # --- 3. Optional MSMC-IM: time-resolved migration rate m(t) and M(t) ---
    if [[ ${RUN_IM} == 1 && ! -s ${OUTPUT_DIR}/msmcim/${PAIR}.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt ]]; then
      mkdir -p ${OUTPUT_DIR}/msmcim
      python3 ${MSMCIM}/MSMC_IM.py -mu ${MU} -p "${P_PAR}" -beta 1e-8,1e-6 \
        -o ${OUTPUT_DIR}/msmcim/${PAIR} --printfittingdetails --plotfittingdetails --xlog ${COMB}
      echo "[$(date)] MSMC-IM done for ${PAIR}"
    fi
  done
done

echo "done running msmc2 for ${TAG}"
date
