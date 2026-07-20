#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   PARTITION=GPUp40 ./run_xgb_chunked_by_node.sh CONFIG_SUBDIR NODE_CPUS NODE_GPUS ARRAY_SPEC
#
# Example:
#   PARTITION=GPUp40 ./run_xgb_chunked_by_node.sh Franzosa_IBD_1 64 1 "1"
#
# Notes:
# - One-task-at-a-time GPU workflow
# - One batch job per selected task
# - Batch job requests the GPU
# - Inner srun step also explicitly requests one GPU for reliable step binding

PARTITION="${PARTITION:-GPUA100}"

CONFIG_SUBDIR="${1:?Need CONFIG_SUBDIR}"
NODE_CPUS="${2:?Need NODE_CPUS}"
NODE_GPUS="${3:?Need NODE_GPUS}"
ARRAY_SPEC="${4:?Need ARRAY_SPEC}"

BASE_DIR="${BASE_DIR:-/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code}"
BASH_DIR="${BASE_DIR}/B_BashScripts"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
LOG_DIR="${BASH_DIR}/logs"
WORKER_SCRIPT="${BASH_DIR}/submit_integrative_array_xgb_single_node.sh"

CPUS_PER_TASK="${CPUS_PER_TASK:-16}"
GPUS_PER_TASK="${GPUS_PER_TASK:-1}"
GPU_GRES_TYPE="${GPU_GRES_TYPE:-gpu}"
USE_GPU="${USE_GPU:-TRUE}"

mkdir -p "$LOG_DIR"

TASKS_BY_CPU=$(( NODE_CPUS / CPUS_PER_TASK ))
TASKS_BY_GPU=$(( NODE_GPUS / GPUS_PER_TASK ))

if (( TASKS_BY_CPU < 1 )); then
    echo "ERROR: NODE_CPUS ($NODE_CPUS) is smaller than CPUS_PER_TASK ($CPUS_PER_TASK)"
    exit 1
fi

if (( TASKS_BY_GPU < 1 )); then
    echo "ERROR: NODE_GPUS ($NODE_GPUS) is smaller than GPUS_PER_TASK ($GPUS_PER_TASK)"
    exit 1
fi

# Force one task per submitted batch job for the current workflow.
TASKS_PER_NODE=1

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

echo "Submitting XGB GPU tasks"
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
echo "NODE_GPUS:       $NODE_GPUS"
echo "CPUS_PER_TASK:   $CPUS_PER_TASK"
echo "GPUS_PER_TASK:   $GPUS_PER_TASK"
echo "TASKS_PER_NODE:  $TASKS_PER_NODE"
echo "GPU_GRES_TYPE:   $GPU_GRES_TYPE"
echo "USE_GPU:         $USE_GPU"
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

    total_gpus_this_chunk=$(( chunk_size * GPUS_PER_TASK ))

    echo "Submitting chunk $chunk_num: $chunk_size task(s)"
    echo "  TASK_IDS = ${chunk[*]}"
    echo "  GPUs     = ${total_gpus_this_chunk}"

    sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --job-name=xgbgpu_${CONFIG_SUBDIR}_c${chunk_num}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --gres=${GPU_GRES_TYPE}:${total_gpus_this_chunk}
#SBATCH --output=${LOG_DIR}/xgbgpu_${CONFIG_SUBDIR}_chunk${chunk_num}_%j.out
#SBATCH --error=${LOG_DIR}/xgbgpu_${CONFIG_SUBDIR}_chunk${chunk_num}_%j.err

set -euo pipefail

BASE_DIR="${BASE_DIR}"
WORKER_SCRIPT="${WORKER_SCRIPT}"
CONFIG_SUBDIR="${CONFIG_SUBDIR}"
LOCAL_CPUS_PER_TASK="${CPUS_PER_TASK}"
LOCAL_GPUS_PER_TASK="${GPUS_PER_TASK}"
USE_GPU="${USE_GPU}"
GPU_GRES_TYPE="${GPU_GRES_TYPE}"
TASK_IDS=( ${task_ids_str} )

echo "SLURM_JOB_ID=\$SLURM_JOB_ID"
echo "SLURM_NODELIST=\$SLURM_NODELIST"
echo "CONFIG_SUBDIR=\$CONFIG_SUBDIR"
echo "WORKER_SCRIPT=\$WORKER_SCRIPT"
echo "Running ${chunk_size} GPU task(s) on one node"
echo "Task IDs: \${TASK_IDS[*]}"
echo "CUDA_VISIBLE_DEVICES(before srun)=\${CUDA_VISIBLE_DEVICES:-unset}"
echo "SLURM_JOB_GPUS(before srun)=\${SLURM_JOB_GPUS:-unset}"
nvidia-smi -L || true
echo

if (( \${#TASK_IDS[@]} != 1 )); then
    echo "ERROR: This launcher expects exactly 1 task per batch job."
    echo "Got TASK_IDS: \${TASK_IDS[*]}"
    exit 1
fi

tid="\${TASK_IDS[0]}"
echo "[LAUNCH] task=\$tid gpu-step via srun"

echo "[LAUNCH] task=\$tid direct worker launch"

export CONFIG_SUBDIR="\$CONFIG_SUBDIR"
export BASE_DIR="\$BASE_DIR"
export SLURM_ARRAY_TASK_ID="\$tid"
export LOCAL_CPUS_PER_TASK="\$LOCAL_CPUS_PER_TASK"
export LOCAL_GPUS_PER_TASK="\$LOCAL_GPUS_PER_TASK"
export LOCAL_JOB_ID="\${SLURM_JOB_ID}_\$tid"
export LOCAL_NODE_NAME="\$(hostname)"
export LOCAL_MEM="node_alloc"
export USE_GPU="\$USE_GPU"
export GPU_DEVICE="cuda"

bash "\$WORKER_SCRIPT"
EOF
done

echo
echo "Done."
echo "Submitted $submitted_total selected task(s) across $chunk_num job(s)."
