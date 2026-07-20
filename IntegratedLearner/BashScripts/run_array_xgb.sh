#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_array_xgb.sh <CONFIG_SUBDIR> [MAX_CONCURRENT] [ARRAY_SPEC]
# Examples:
#   ./run_array_xgb.sh Erawijantari 10
#   ./run_array_xgb.sh Erawijantari 10 0-31

CONFIG_SUBDIR="${1:?ERROR: provide CONFIG_SUBDIR (e.g., Erawijantari)}"
MAX_CONCURRENT="${2:-10}"
ARRAY_SPEC="${3:-}"

BASE_DIR="/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"

if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "ERROR: CONFIG_DIR does not exist: $CONFIG_DIR"
  exit 1
fi

N=$(ls -1 "$CONFIG_DIR"/*.txt 2>/dev/null | wc -l | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No .txt config files found in: $CONFIG_DIR"
  exit 1
fi

# xgb-only: metab(2) * taxa(2) * reduction(2) * models(1) = 8
PER_CONFIG=8
TOTAL=$((N * PER_CONFIG))

mkdir -p "${BASE_DIR}/B_BashScripts/logs"

if [[ -n "$ARRAY_SPEC" ]]; then
  ARRAY_ARG="--array=${ARRAY_SPEC}%${MAX_CONCURRENT}"
else
  ARRAY_ARG="--array=0-$((TOTAL-1))%${MAX_CONCURRENT}"
fi

echo "Submitting XGB array for CONFIG_SUBDIR=${CONFIG_SUBDIR}"
echo "Configs found: ${N}"
echo "PER_CONFIG: ${PER_CONFIG}"
echo "TOTAL_TASKS: ${TOTAL}"
echo "MAX_CONCURRENT: ${MAX_CONCURRENT}"
echo "ARRAY: ${ARRAY_ARG}"

sbatch --export=ALL,CONFIG_SUBDIR="${CONFIG_SUBDIR}",BASE_DIR="${BASE_DIR}" \
  ${ARRAY_ARG} \
  submit_integrative_array_xgb.sh
