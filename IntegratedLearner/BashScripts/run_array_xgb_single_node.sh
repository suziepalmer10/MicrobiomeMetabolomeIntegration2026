#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_array_xgb_single_node.sh <CONFIG_SUBDIR> [MAX_CONCURRENT] [ARRAY_SPEC] [CPUS_PER_TASK]
#
# Examples:
#   ./run_array_xgb_single_node.sh Erawijantari 10
#   ./run_array_xgb_single_node.sh Erawijantari 10 0-31
#   ./run_array_xgb_single_node.sh Erawijantari 4 0-7,12,20-23 16
#
# Notes:
# - Runs the entire former Slurm array locally on one node.
# - MAX_CONCURRENT = number of simultaneous array tasks to run.
# - CPUS_PER_TASK = threads assigned to each task.

CONFIG_SUBDIR="${1:?ERROR: provide CONFIG_SUBDIR (e.g., Erawijantari)}"
MAX_CONCURRENT="${2:-10}"
ARRAY_SPEC="${3:-}"
CPUS_PER_TASK="${4:-16}"

BASE_DIR="/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
WORKER_SCRIPT="${BASE_DIR}/B_BashScripts/submit_integrative_array_xgb_single_node.sh"
LOG_DIR="${BASE_DIR}/B_BashScripts/logs"

if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "ERROR: CONFIG_DIR does not exist: $CONFIG_DIR"
  exit 1
fi

if [[ ! -f "$WORKER_SCRIPT" ]]; then
  echo "ERROR: Worker script not found: $WORKER_SCRIPT"
  exit 1
fi

N=$(find "$CONFIG_DIR" -maxdepth 1 -type f -name '*.txt' | wc -l | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No .txt config files found in: $CONFIG_DIR"
  exit 1
fi

# xgb-only: metab(2) * taxa(2) * reduction(2) * models(1) = 8
PER_CONFIG=8
TOTAL=$((N * PER_CONFIG))

mkdir -p "$LOG_DIR"

expand_array_spec() {
  local spec="$1"
  local total="$2"
  local part start end i
  local -a out=()

  if [[ -z "$spec" ]]; then
    for ((i=0; i<total; i++)); do
      out+=("$i")
    done
  else
    IFS=',' read -r -a parts <<< "$spec"
    for part in "${parts[@]}"; do
      if [[ "$part" =~ ^[0-9]+$ ]]; then
        if (( part < 0 || part >= total )); then
          echo "ERROR: task index out of range: $part (valid: 0-$((total-1)))" >&2
          exit 1
        fi
        out+=("$part")
      elif [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        start="${BASH_REMATCH[1]}"
        end="${BASH_REMATCH[2]}"
        if (( start > end )); then
          echo "ERROR: invalid range: $part" >&2
          exit 1
        fi
        if (( start < 0 || end >= total )); then
          echo "ERROR: range out of bounds: $part (valid: 0-$((total-1)))" >&2
          exit 1
        fi
        for ((i=start; i<=end; i++)); do
          out+=("$i")
        done
      else
        echo "ERROR: unsupported ARRAY_SPEC token: $part" >&2
        echo "Use forms like: 0-31 or 0-7,12,20-23" >&2
        exit 1
      fi
    done
  fi

  printf '%s\n' "${out[@]}"
}

wait_for_slot() {
  local max_jobs="$1"
  while true; do
    local running
    running=$(jobs -rp | wc -l | tr -d ' ')
    if (( running < max_jobs )); then
      break
    fi
    sleep 1
  done
}

TASK_IDS=()
while IFS= read -r tid; do
  TASK_IDS+=("$tid")
done < <(expand_array_spec "$ARRAY_SPEC" "$TOTAL")

if [[ "${#TASK_IDS[@]}" -eq 0 ]]; then
  echo "ERROR: No tasks selected."
  exit 1
fi

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
LOCAL_JOB_ID="xgb_single_node_${CONFIG_SUBDIR}_${RUN_STAMP}"
LOCAL_NODE_NAME="$(hostname)"

echo "Running XGB array on one node"
echo "CONFIG_SUBDIR:   ${CONFIG_SUBDIR}"
echo "Configs found:   ${N}"
echo "PER_CONFIG:      ${PER_CONFIG}"
echo "TOTAL_TASKS:     ${TOTAL}"
echo "SELECTED_TASKS:  ${#TASK_IDS[@]}"
echo "MAX_CONCURRENT:  ${MAX_CONCURRENT}"
echo "CPUS_PER_TASK:   ${CPUS_PER_TASK}"
echo "NODE:            ${LOCAL_NODE_NAME}"
echo "JOB_ID:          ${LOCAL_JOB_ID}"
if [[ -n "$ARRAY_SPEC" ]]; then
  echo "ARRAY_SPEC:      ${ARRAY_SPEC}"
else
  echo "ARRAY_SPEC:      full range 0-$((TOTAL-1))"
fi
echo

pids=()
pid_to_task=()

for task_id in "${TASK_IDS[@]}"; do
  wait_for_slot "$MAX_CONCURRENT"

  stdout_log="${LOG_DIR}/${LOCAL_JOB_ID}_task_${task_id}.out"
  stderr_log="${LOG_DIR}/${LOCAL_JOB_ID}_task_${task_id}.err"

  (
    set -euo pipefail
    export CONFIG_SUBDIR="$CONFIG_SUBDIR"
    export BASE_DIR="$BASE_DIR"
    export SLURM_ARRAY_TASK_ID="$task_id"
    export SLURM_ARRAY_TASK_MIN=0
    export SLURM_ARRAY_TASK_MAX=$((TOTAL - 1))
    export SLURM_ARRAY_TASK_COUNT="$TOTAL"

    export LOCAL_CPUS_PER_TASK="$CPUS_PER_TASK"
    export LOCAL_JOB_ID="$LOCAL_JOB_ID"
    export LOCAL_NODE_NAME="$LOCAL_NODE_NAME"
    export LOCAL_MEM="single_node"

    bash "$WORKER_SCRIPT"
  ) >"$stdout_log" 2>"$stderr_log" &

  pid=$!
  pids+=("$pid")
  pid_to_task+=("${pid}:${task_id}")

  echo "Launched task ${task_id}"
  echo "  stdout: $stdout_log"
  echo "  stderr: $stderr_log"
done

echo
echo "Waiting for all tasks to finish..."

fail=0

for entry in "${pid_to_task[@]}"; do
  pid="${entry%%:*}"
  task_id="${entry##*:}"
  if wait "$pid"; then
    echo "Task ${task_id} completed successfully"
  else
    echo "Task ${task_id} FAILED"
    fail=1
  fi
done

echo
if (( fail != 0 )); then
  echo "One or more tasks failed."
  echo "Check logs in: $LOG_DIR"
  exit 1
fi

echo "All selected tasks completed successfully."
echo "Logs are in: $LOG_DIR"
