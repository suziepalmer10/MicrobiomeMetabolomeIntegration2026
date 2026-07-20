#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   PARTITION=CPU ./run_chunked_xgb_multi_cpu.sh CONFIG_SUBDIR NODE_CPUS ARRAY_SPEC
#
# Example:
#   PARTITION=CPU ./run_chunked_xgb_multi_cpu.sh Wang 64 "1,3,5,7"

PARTITION="${PARTITION:-CPU}"

CONFIG_SUBDIR="${1:?Need CONFIG_SUBDIR}"
NODE_CPUS="${2:?Need NODE_CPUS}"
ARRAY_SPEC="${3:?Need ARRAY_SPEC}"

BASE_DIR="${BASE_DIR:-/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code}"
BASH_DIR="${BASE_DIR}/B_BashScripts"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
LOG_DIR="${BASH_DIR}/logs"
WORKER_SCRIPT="${BASH_DIR}/submit_integrative_array_xgb_single_node_cpu.sh"

CPUS_PER_TASK="${CPUS_PER_TASK:-16}"
TASKS_PER_NODE=1
USE_GPU="${USE_GPU:-FALSE}"
GPU_DEVICE="${GPU_DEVICE:-cpu}"
LOCAL_GPUS_PER_TASK=0

mkdir -p "$LOG_DIR"

TASKS_BY_CPU=$(( NODE_CPUS / CPUS_PER_TASK ))

if (( TASKS_BY_CPU < 1 )); then
    echo "ERROR: NODE_CPUS ($NODE_CPUS) is smaller than CPUS_PER_TASK ($CPUS_PER_TASK)"
    exit 1
fi

if [[ ! -d "$CONFIG_DIR" ]]; then
    echo "ERROR: Config directory not found: $CONFIG_DIR"
    exit 1
fi

if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Worker script not found: $WORKER_SCRIPT"
    exit 1
fi

models=(xgb)
metab_transforms=(log2 none)
taxa_transforms=(clr none)
reduction_pairs=(limma none)

PER_CONFIG=$(( ${#models[@]} * ${#metab_transforms[@]} * ${#taxa_transforms[@]} * ${#reduction_pairs[@]} ))

mapfile -t CONFIG_FILES < <(find "$CONFIG_DIR" -maxdepth 1 -type f -name '*.txt' | sort)
NUM_CONFIGS="${#CONFIG_FILES[@]}"

if (( NUM_CONFIGS == 0 )); then
    echo "ERROR: No config files found in $CONFIG_DIR"
    exit 1
fi

TOTAL_TASKS=$(( NUM_CONFIGS * PER_CONFIG ))

echo "Submitting XGB CPU tasks"
echo "BASE_DIR:        $BASE_DIR"
echo "BASH_DIR:        $BASH_DIR"
echo "WORKER_SCRIPT:   $WORKER_SCRIPT"
echo "CONFIG_SUBDIR:   $CONFIG_SUBDIR"
echo "CONFIG_DIR:      $CONFIG_DIR"
echo "Configs found:   $NUM_CONFIGS"
echo "PER_CONFIG:      $PER_CONFIG"
echo "TOTAL_TASKS:     $TOTAL_TASKS"
echo "ARRAY_SPEC:      $ARRAY_SPEC"
echo "PARTITION:       $PARTITION"
echo "NODE_CPUS:       $NODE_CPUS"
echo "CPUS_PER_TASK:   $CPUS_PER_TASK"
echo "TASKS_PER_NODE:  $TASKS_PER_NODE"
echo "USE_GPU:         $USE_GPU"
echo "GPU_DEVICE:      $GPU_DEVICE"
echo

expand_array_spec() {
    local spec="$1"
    local IFS=','
    local parts=()
    read -ra parts <<< "$spec"

    for part in "${parts[@]}"; do
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            local start="${BASH_REMATCH[1]}"
            local end="${BASH_REMATCH[2]}"
            local i
            for (( i=start; i<=end; i++ )); do
                echo "$i"
            done
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            echo "$part"
        else
            echo "ERROR: Invalid ARRAY_SPEC token: $part" >&2
            exit 1
        fi
    done
}

mapfile -t SELECTED_TASKS < <(expand_array_spec "$ARRAY_SPEC")

VALIDATED_TASKS=()
for tid in "${SELECTED_TASKS[@]}"; do
    if (( tid < 0 || tid >= TOTAL_TASKS )); then
        echo "WARNING: Skipping out-of-range task id $tid (valid range: 0-$((TOTAL_TASKS - 1)))"
        continue
    fi
    VALIDATED_TASKS+=("$tid")
done
SELECTED_TASKS=("${VALIDATED_TASKS[@]}")

NUM_SELECTED="${#SELECTED_TASKS[@]}"
echo "SELECTED_TASKS:  $NUM_SELECTED"
echo

if (( NUM_SELECTED == 0 )); then
    echo "Nothing to submit."
    exit 0
fi

chunk_num=0
submitted_total=0

for (( offset=0; offset<NUM_SELECTED; offset+=TASKS_PER_NODE )); do
    chunk_num=$((chunk_num + 1))
    chunk=( "${SELECTED_TASKS[@]:offset:TASKS_PER_NODE}" )
    chunk_size="${#chunk[@]}"
    submitted_total=$((submitted_total + chunk_size))
    task_ids_str="${chunk[*]}"

    echo "Submitting chunk $chunk_num: $chunk_size task(s)"
    echo "  TASK_IDS = ${chunk[*]}"
    echo "  CPUs     = ${CPUS_PER_TASK}"

    sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --job-name=xgbcpu_${CONFIG_SUBDIR}_c${chunk_num}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${LOG_DIR}/xgbcpu_${CONFIG_SUBDIR}_chunk${chunk_num}_%j.out
#SBATCH --error=${LOG_DIR}/xgbcpu_${CONFIG_SUBDIR}_chunk${chunk_num}_%j.err

set -euo pipefail

BASE_DIR="${BASE_DIR}"
WORKER_SCRIPT="${WORKER_SCRIPT}"
CONFIG_SUBDIR="${CONFIG_SUBDIR}"
LOCAL_CPUS_PER_TASK="${CPUS_PER_TASK}"
LOCAL_GPUS_PER_TASK="${LOCAL_GPUS_PER_TASK}"
USE_GPU="${USE_GPU}"
GPU_DEVICE="${GPU_DEVICE}"
TASK_IDS=( ${task_ids_str} )

echo "SLURM_JOB_ID=\$SLURM_JOB_ID"
echo "SLURM_NODELIST=\$SLURM_NODELIST"
echo "CONFIG_SUBDIR=\$CONFIG_SUBDIR"
echo "WORKER_SCRIPT=\$WORKER_SCRIPT"
echo "Running ${chunk_size} CPU task(s) on one node"
echo "Task IDs: \${TASK_IDS[*]}"
echo

if (( \${#TASK_IDS[@]} != 1 )); then
    echo "ERROR: This launcher expects exactly 1 task per batch job."
    echo "Got TASK_IDS: \${TASK_IDS[*]}"
    exit 1
fi

tid="\${TASK_IDS[0]}"
echo "[LAUNCH] task=\$tid cpu worker launch"

export CONFIG_SUBDIR="\$CONFIG_SUBDIR"
export BASE_DIR="\$BASE_DIR"
export SLURM_ARRAY_TASK_ID="\$tid"
export LOCAL_CPUS_PER_TASK="\$LOCAL_CPUS_PER_TASK"
export LOCAL_GPUS_PER_TASK="\$LOCAL_GPUS_PER_TASK"
export LOCAL_JOB_ID="\${SLURM_JOB_ID}_\$tid"
export LOCAL_NODE_NAME="\$(hostname)"
export LOCAL_MEM="node_alloc"
export USE_GPU="\$USE_GPU"
export GPU_DEVICE="\$GPU_DEVICE"

bash "\$WORKER_SCRIPT"
EOF
done

echo
echo "Done."
echo "Submitted $submitted_total selected task(s) across $chunk_num job(s)."
