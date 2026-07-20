#!/bin/bash
set -euo pipefail

# Usage:
#   ./run_small_chunked_by_node.sh <CONFIG_SUBDIR> [TASKS_PER_NODE] [ARRAY_SPEC]
#
# Examples:
#   ./run_small_chunked_by_node.sh Wang 64
#   ./run_small_chunked_by_node.sh Wang 64 "1,3,5,7"
#   ./run_small_chunked_by_node.sh Wang 64 "0-15"

CONFIG_SUBDIR="${1:?ERROR: provide CONFIG_SUBDIR}"
TASKS_PER_NODE="${2:-64}"
ARRAY_SPEC="${3:-}"

BASE_DIR="/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
LOG_DIR="${BASE_DIR}/B_BashScripts/logs"

if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "ERROR: CONFIG_DIR does not exist: $CONFIG_DIR"
  exit 1
fi

N=$(find "$CONFIG_DIR" -maxdepth 1 -type f -name '*.txt' | wc -l | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No .txt config files found in: $CONFIG_DIR"
  exit 1
fi

# small models = enet + rf
# models(2) * metab(2) * taxa(2) * reduction(2) = 16
PER_CONFIG=16
TOTAL=$((N * PER_CONFIG))

mkdir -p "$LOG_DIR"

expand_array_spec() {
  local spec="$1"
  local out=()
  local part start end i
  local parts=()

  IFS=',' read -ra parts <<< "$spec"
  for part in "${parts[@]}"; do
    if [[ "$part" =~ ^[0-9]+-[0-9]+$ ]]; then
      start="${part%-*}"
      end="${part#*-}"
      if (( start > end )); then
        echo "ERROR: Invalid range in ARRAY_SPEC: $part" >&2
        exit 1
      fi
      for ((i=start; i<=end; i++)); do
        out+=("$i")
      done
    elif [[ "$part" =~ ^[0-9]+$ ]]; then
      out+=("$part")
    else
      echo "ERROR: Invalid ARRAY_SPEC component: $part" >&2
      exit 1
    fi
  done

  printf '%s\n' "${out[@]}" | awk '!seen[$0]++'
}

if [[ -n "$ARRAY_SPEC" ]]; then
  mapfile -t TASK_IDS < <(expand_array_spec "$ARRAY_SPEC")
else
  mapfile -t TASK_IDS < <(seq 0 $((TOTAL - 1)))
fi

if [[ "${#TASK_IDS[@]}" -eq 0 ]]; then
  echo "ERROR: No task IDs resolved."
  exit 1
fi

for tid in "${TASK_IDS[@]}"; do
  if (( tid < 0 || tid >= TOTAL )); then
    echo "ERROR: Task ID $tid out of range 0..$((TOTAL-1))"
    exit 1
  fi
done

echo "Submitting SMALL node-chunked jobs for CONFIG_SUBDIR=${CONFIG_SUBDIR}"
echo "Configs found: ${N}"
echo "PER_CONFIG: ${PER_CONFIG}"
echo "TOTAL_TASKS: ${TOTAL}"
echo "TASKS_PER_NODE: ${TASKS_PER_NODE}"
echo "Requested task IDs: ${TASK_IDS[*]}"

chunk=()
chunk_idx=0

submit_chunk() {
  local ids_display ids_export n_tasks
  ids_display="$(IFS=,; echo "${chunk[*]}")"
  ids_export="$(printf '%s ' "${chunk[@]}")"
  ids_export="${ids_export% }"
  n_tasks="${#chunk[@]}"

  echo
  echo "Submitting chunk ${chunk_idx} with task IDs: ${ids_display}"

  sbatch \
    --ntasks="${n_tasks}" \
    --export=ALL,CONFIG_SUBDIR="${CONFIG_SUBDIR}",BASE_DIR="${BASE_DIR}",TASK_IDS="${ids_export}" \
    "${BASE_DIR}/B_BashScripts/submit_integrative_array_small.sh"
}

for tid in "${TASK_IDS[@]}"; do
  chunk+=("$tid")
  if (( ${#chunk[@]} == TASKS_PER_NODE )); then
    submit_chunk
    chunk=()
    chunk_idx=$((chunk_idx + 1))
  fi
done

if (( ${#chunk[@]} > 0 )); then
  submit_chunk
fi
