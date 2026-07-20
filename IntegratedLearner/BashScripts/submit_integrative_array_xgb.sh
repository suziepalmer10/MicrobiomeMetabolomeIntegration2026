#!/bin/bash
#SBATCH --partition=super
#SBATCH --time=99:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=480G
#SBATCH --hint=nomultithread
#SBATCH --mail-type=ALL
#SBATCH --mail-user=suzette.palmer@utsouthwestern.edu
#SBATCH --output=/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code/B_BashScripts/logs/%x_%j.out
#SBATCH --error=/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code/B_BashScripts/logs/%x_%j.err

set -euo pipefail

: "${CONFIG_SUBDIR:?ERROR: CONFIG_SUBDIR not set}"
: "${ARRAY_SPEC:?ERROR: ARRAY_SPEC not set (example: 0-3 or 1,3,5,7)}"

BASE_DIR="${BASE_DIR:-/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code}"
PIPELINE_DIR="${BASE_DIR}/A_Integrative_Pipeline_Scripts"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
LOG_DIR="${BASE_DIR}/B_BashScripts/logs"

TASKS_PER_NODE="${TASKS_PER_NODE:-4}"
CPUS_PER_TASK_LOCAL="${CPUS_PER_TASK_LOCAL:-16}"
MEM_PER_TASK_LOCAL="${MEM_PER_TASK_LOCAL:-120G}"

if [[ ! -d "$PIPELINE_DIR" ]]; then
  echo "ERROR: PIPELINE_DIR not found: $PIPELINE_DIR"
  exit 1
fi

if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "ERROR: CONFIG_DIR not found: $CONFIG_DIR"
  exit 1
fi

mkdir -p "$LOG_DIR"
cd "$PIPELINE_DIR"

models=(xgb)
metab_transforms=(log2 none)
taxa_transforms=(clr none)
reduction_pairs=(limma none)

PER_CONFIG=$(( ${#models[@]} * ${#metab_transforms[@]} * ${#taxa_transforms[@]} * ${#reduction_pairs[@]} ))

mapfile -t configs < <(find "$CONFIG_DIR" -maxdepth 1 -type f -name '*.txt' | sort)
NUM_CONFIGS=${#configs[@]}

if [[ "$NUM_CONFIGS" -eq 0 ]]; then
  echo "ERROR: No config files found in $CONFIG_DIR"
  exit 1
fi

TOTAL_TASKS=$(( NUM_CONFIGS * PER_CONFIG ))

expand_task_spec() {
  local spec="$1"
  local part start end i
  local -a expanded=()

  IFS=',' read -r -a parts <<< "$spec"
  for part in "${parts[@]}"; do
    part="${part//[[:space:]]/}"
    [[ -z "$part" ]] && continue

    if [[ "$part" =~ ^[0-9]+-[0-9]+$ ]]; then
      start="${part%-*}"
      end="${part#*-}"
      if (( start > end )); then
        echo "ERROR: Bad range in ARRAY_SPEC: $part" >&2
        exit 1
      fi
      for (( i=start; i<=end; i++ )); do
        expanded+=("$i")
      done
    elif [[ "$part" =~ ^[0-9]+$ ]]; then
      expanded+=("$part")
    else
      echo "ERROR: Invalid token in ARRAY_SPEC: $part" >&2
      exit 1
    fi
  done

  printf '%s\n' "${expanded[@]}"
}

mapfile -t requested_task_ids < <(expand_task_spec "$ARRAY_SPEC")

if [[ "${#requested_task_ids[@]}" -eq 0 ]]; then
  echo "ERROR: ARRAY_SPEC expanded to zero tasks"
  exit 1
fi

if (( ${#requested_task_ids[@]} > TASKS_PER_NODE )); then
  echo "ERROR: ARRAY_SPEC expanded to ${#requested_task_ids[@]} tasks, but TASKS_PER_NODE=$TASKS_PER_NODE"
  echo "Submit at most $TASKS_PER_NODE task IDs per node-chunk job."
  exit 1
fi

declare -A seen=()
valid_task_ids=()
for task_id in "${requested_task_ids[@]}"; do
  if (( task_id < 0 || task_id >= TOTAL_TASKS )); then
    echo "ERROR: TASK_ID=$task_id out of range 0..$((TOTAL_TASKS-1))"
    exit 1
  fi
  if [[ -z "${seen[$task_id]:-}" ]]; then
    valid_task_ids+=("$task_id")
    seen["$task_id"]=1
  fi
done

echo "==== NODE-CHUNK JOB INFO ===="
echo "JobID:           ${SLURM_JOB_ID:-unknown}"
echo "Node(s):         ${SLURM_JOB_NODELIST:-unknown}"
echo "CONFIG_SUBDIR:   $CONFIG_SUBDIR"
echo "BASE_DIR:        $BASE_DIR"
echo "ARRAY_SPEC:      $ARRAY_SPEC"
echo "Expanded tasks:  ${valid_task_ids[*]}"
echo "TOTAL_TASKS:     $TOTAL_TASKS"
echo "TASKS_PER_NODE:  $TASKS_PER_NODE"
echo "SLURM_NTASKS:    ${SLURM_NTASKS:-unknown}"
echo "CPUS_PER_TASK:   ${SLURM_CPUS_PER_TASK:-unknown}"
echo "============================="

pids=()
task_logs=()

for task_id in "${valid_task_ids[@]}"; do
  task_log="${LOG_DIR}/xgb_${CONFIG_SUBDIR}_job${SLURM_JOB_ID:-nojid}_task${task_id}.out"
  task_err="${LOG_DIR}/xgb_${CONFIG_SUBDIR}_job${SLURM_JOB_ID:-nojid}_task${task_id}.err"

  echo "Launching TASK_ID=$task_id"
  echo "  stdout -> $task_log"
  echo "  stderr -> $task_err"

  srun --exclusive -N1 -n1 -c "$CPUS_PER_TASK_LOCAL" --cpu-bind=cores \
    env \
      TASK_ID="$task_id" \
      CONFIG_SUBDIR="$CONFIG_SUBDIR" \
      BASE_DIR="$BASE_DIR" \
      LOCAL_CPUS_PER_TASK="$CPUS_PER_TASK_LOCAL" \
      LOCAL_JOB_ID="${SLURM_JOB_ID:-single_node}" \
      LOCAL_NODE_NAME="${SLURM_JOB_NODELIST:-$(hostname)}" \
      LOCAL_MEM="$MEM_PER_TASK_LOCAL" \
    bash submit_integrative_array_xgb_single_node.sh \
    >"$task_log" 2>"$task_err" &

  pids+=("$!")
  task_logs+=("$task_id:$task_log:$task_err")
done

rc=0

for idx in "${!pids[@]}"; do
  pid="${pids[$idx]}"
  meta="${task_logs[$idx]}"
  task_id="${meta%%:*}"
  rest="${meta#*:}"
  task_log="${rest%%:*}"
  task_err="${rest#*:}"

  if wait "$pid"; then
    echo "TASK_ID=$task_id finished successfully"
  else
    echo "TASK_ID=$task_id failed"
    echo "  stdout: $task_log"
    echo "  stderr: $task_err"
    rc=1
  fi
done

exit "$rc"
