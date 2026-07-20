#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   PARTITION=GPU4v100 ./run_xgb_pack_on_one_node.sh CONFIG_SUBDIR NODE_CPUS NODE_GPUS ARRAY_SPEC
#
# Example:
#   PARTITION=GPU4v100 ./run_xgb_pack_on_one_node.sh Franzosa_CD_2 64 4 "0-7"
#   PARTITION=GPU4v100 ./run_xgb_pack_on_one_node.sh Franzosa_IBD_2 64 4 "0,1,2,3"

PARTITION="${PARTITION:-GPU4v100}"

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

if [[ ! -d "$CONFIG_DIR" ]]; then
    echo "ERROR: Config directory not found: $CONFIG_DIR"
    exit 1
fi

if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Worker script not found: $WORKER_SCRIPT"
    exit 1
fi

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

MAX_CONCURRENT_TASKS=$(( TASKS_BY_CPU < TASKS_BY_GPU ? TASKS_BY_CPU : TASKS_BY_GPU ))

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

if (( NUM_SELECTED == 0 )); then
    echo "Nothing to run."
    exit 0
fi

echo "Submitting packed XGB GPU job"
echo "CONFIG_SUBDIR:         $CONFIG_SUBDIR"
echo "CONFIG_DIR:            $CONFIG_DIR"
echo "TOTAL_TASKS:           $TOTAL_TASKS"
echo "ARRAY_SPEC:            $ARRAY_SPEC"
echo "SELECTED_TASKS:        ${SELECTED_TASKS[*]}"
echo "PARTITION:             $PARTITION"
echo "NODE_CPUS:             $NODE_CPUS"
echo "NODE_GPUS:             $NODE_GPUS"
echo "CPUS_PER_TASK:         $CPUS_PER_TASK"
echo "GPUS_PER_TASK:         $GPUS_PER_TASK"
echo "MAX_CONCURRENT_TASKS:  $MAX_CONCURRENT_TASKS"
echo

if (( NUM_SELECTED > MAX_CONCURRENT_TASKS )); then
    echo "ERROR: Selected $NUM_SELECTED tasks, but node can only run $MAX_CONCURRENT_TASKS concurrently."
    echo "Either reduce ARRAY_SPEC or use a bigger node / smaller per-task resources."
    exit 1
fi

task_ids_str="${SELECTED_TASKS[*]}"

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --job-name=xgbpack_${CONFIG_SUBDIR}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${NODE_CPUS}
#SBATCH --gres=${GPU_GRES_TYPE}:${NUM_SELECTED}
#SBATCH --output=${LOG_DIR}/xgbpack_${CONFIG_SUBDIR}_%j.out
#SBATCH --error=${LOG_DIR}/xgbpack_${CONFIG_SUBDIR}_%j.err

set -euo pipefail

BASE_DIR="${BASE_DIR}"
WORKER_SCRIPT="${WORKER_SCRIPT}"
CONFIG_SUBDIR="${CONFIG_SUBDIR}"
LOCAL_CPUS_PER_TASK="${CPUS_PER_TASK}"
LOCAL_GPUS_PER_TASK="${GPUS_PER_TASK}"
USE_GPU="${USE_GPU}"
TASK_IDS=( ${task_ids_str} )

echo "SLURM_JOB_ID=\$SLURM_JOB_ID"
echo "SLURM_NODELIST=\$SLURM_NODELIST"
echo "CONFIG_SUBDIR=\$CONFIG_SUBDIR"
echo "TASK_IDS=\${TASK_IDS[*]}"
echo "CUDA_VISIBLE_DEVICES(before launch)=\${CUDA_VISIBLE_DEVICES:-unset}"
echo "SLURM_JOB_GPUS=\${SLURM_JOB_GPUS:-unset}"
nvidia-smi -L || true
echo

pids=()
gpu_idx=0

for tid in "\${TASK_IDS[@]}"; do
    echo "[LAUNCH] task=\$tid gpu=\$gpu_idx"

    (
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
        export CUDA_VISIBLE_DEVICES="\$gpu_idx"
        bash "\$WORKER_SCRIPT"
    ) &
    pids+=(\$!)
    gpu_idx=\$((gpu_idx + 1))
done

rc=0
for pid in "\${pids[@]}"; do
    if ! wait "\$pid"; then
        rc=1
    fi
done

exit \$rc
EOF
