#!/bin/bash
set -euo pipefail

: "${CONFIG_SUBDIR:?ERROR: CONFIG_SUBDIR not set}"
BASE_DIR="${BASE_DIR:-/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code}"

if [[ -n "${TASK_ID:-}" ]]; then
  RUN_TASK_ID="$TASK_ID"
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  RUN_TASK_ID="$SLURM_ARRAY_TASK_ID"
else
  echo "ERROR: neither TASK_ID nor SLURM_ARRAY_TASK_ID is set"
  exit 1
fi

LOCAL_CPUS_PER_TASK="${LOCAL_CPUS_PER_TASK:-${SLURM_CPUS_PER_TASK:-16}}"
LOCAL_GPUS_PER_TASK="${LOCAL_GPUS_PER_TASK:-1}"
LOCAL_JOB_ID="${LOCAL_JOB_ID:-${SLURM_JOB_ID:-single_node}}"
LOCAL_NODE_NAME="${LOCAL_NODE_NAME:-${SLURMD_NODENAME:-$(hostname)}}"
LOCAL_MEM="${LOCAL_MEM:-unknown}"

USE_GPU="${USE_GPU:-FALSE}"
GPU_DEVICE="${GPU_DEVICE:-cuda}"

PIPELINE_DIR="${BASE_DIR}/A_Integrative_Pipeline_Scripts"
CONFIG_DIR="${BASE_DIR}/C_Configuration_Files/${CONFIG_SUBDIR}"
LOG_DIR="${BASE_DIR}/B_BashScripts/logs"

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

export OMP_NUM_THREADS="$LOCAL_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$LOCAL_CPUS_PER_TASK"
export MKL_NUM_THREADS="$LOCAL_CPUS_PER_TASK"
export VECLIB_MAXIMUM_THREADS="$LOCAL_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$LOCAL_CPUS_PER_TASK"
export R_PARALLEL_NUM_WORKERS="$LOCAL_CPUS_PER_TASK"

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

if (( RUN_TASK_ID < 0 || RUN_TASK_ID >= TOTAL_TASKS )); then
  echo "ERROR: TASK_ID=$RUN_TASK_ID out of range 0..$((TOTAL_TASKS-1))"
  exit 1
fi

config_idx=$(( RUN_TASK_ID / PER_CONFIG ))
combo_idx=$(( RUN_TASK_ID % PER_CONFIG ))
config_file="${configs[$config_idx]}"

n_models=${#models[@]}
n_metab=${#metab_transforms[@]}
n_taxa=${#taxa_transforms[@]}
n_red=${#reduction_pairs[@]}

tmp=$combo_idx
red_idx=$(( tmp % n_red ));     tmp=$(( tmp / n_red ))
taxa_idx=$(( tmp % n_taxa ));   tmp=$(( tmp / n_taxa ))
metab_idx=$(( tmp % n_metab )); tmp=$(( tmp / n_metab ))
model_idx=$(( tmp % n_models ))

MODEL_KEY="${models[$model_idx]}"
METAB_TRANSFORM_KEY="${metab_transforms[$metab_idx]}"
TAXA_TRANSFORM_KEY="${taxa_transforms[$taxa_idx]}"
REDUCTION_PAIR_KEY="${reduction_pairs[$red_idx]}"

case "$MODEL_KEY" in
  enet) model_to_run="${PIPELINE_DIR}/model_functions/enet_function.R" ;;
  rf)   model_to_run="${PIPELINE_DIR}/model_functions/rf_function.R" ;;
  xgb)  model_to_run="${PIPELINE_DIR}/model_functions/xgboost_function.R" ;;
  *) echo "ERROR: bad MODEL_KEY=$MODEL_KEY"; exit 1 ;;
esac

case "$METAB_TRANSFORM_KEY" in
  log2) metab_transform="log2" ;;
  none) metab_transform="none" ;;
  *) echo "ERROR: bad METAB_TRANSFORM_KEY=$METAB_TRANSFORM_KEY"; exit 1 ;;
esac

case "$TAXA_TRANSFORM_KEY" in
  clr) taxa_transform="clr" ;;
  none) taxa_transform="none" ;;
  *) echo "ERROR: bad TAXA_TRANSFORM_KEY=$TAXA_TRANSFORM_KEY"; exit 1 ;;
esac

case "$REDUCTION_PAIR_KEY" in
  limma) metab_reduction="limma"; taxa_reduction="wilcox" ;;
  none)  metab_reduction="none";  taxa_reduction="none" ;;
  *) echo "ERROR: bad REDUCTION_PAIR_KEY=$REDUCTION_PAIR_KEY"; exit 1 ;;
esac

if [[ "$MODEL_KEY" == "xgb" ]]; then
  MODEL_USE_GPU="${USE_GPU}"
  MODEL_GPU_DEVICE="${GPU_DEVICE}"
else
  MODEL_USE_GPU="FALSE"
  MODEL_GPU_DEVICE="cpu"
fi

echo "---- LOCAL TASK INFO ----"
echo "JobID:          ${LOCAL_JOB_ID}"
echo "Node:           ${LOCAL_NODE_NAME}"
echo "CPUs/task:      ${LOCAL_CPUS_PER_TASK}"
echo "GPUs/task:      ${LOCAL_GPUS_PER_TASK}"
echo "Mem:            ${LOCAL_MEM}"
echo "BASE_DIR:       $BASE_DIR"
echo "PIPELINE_DIR:   $PIPELINE_DIR"
echo "CONFIG_SUBDIR:  $CONFIG_SUBDIR"
echo "TASK_ID:        $RUN_TASK_ID / $((TOTAL_TASKS-1))"
echo "CONFIG_FILE:    $config_file"
echo "MODEL_KEY:      $MODEL_KEY"
echo "MODEL_SCRIPT:   $model_to_run"
echo "METAB_TRANSFORM:$METAB_TRANSFORM_KEY"
echo "TAXA_TRANSFORM: $TAXA_TRANSFORM_KEY"
echo "REDUCTION_PAIR: $REDUCTION_PAIR_KEY (metab=$metab_reduction, taxa=$taxa_reduction)"
echo "USE_GPU:        $MODEL_USE_GPU"
echo "GPU_DEVICE:     $MODEL_GPU_DEVICE"
echo "SLURM_GPUS_ON_NODE=${SLURM_GPUS_ON_NODE:-unset}"
echo "-------------------------"

source "$config_file"

file_path=${file_path//\'/}
input_file=${input_file//\'/}
type_of_analysis=${type_of_analysis//\'/}
response_variable=${response_variable//\'/}
stratify_variable=${stratify_variable//\'/}
training_proportion=${training_proportion//\'/}
num_repeats=${num_repeats//\'/}
num_folds=${num_folds//\'/}

file_path="$BASE_DIR"

if [[ "$input_file" != /* ]]; then
  input_file="${PIPELINE_DIR}/${input_file}"
fi

if [[ ! -f "$input_file" ]]; then
  echo "ERROR: Resolved input_file not found: $input_file"
  exit 1
fi

if [[ ! -f "$model_to_run" ]]; then
  echo "ERROR: Resolved model_to_run not found from PIPELINE_DIR: $model_to_run"
  exit 1
fi

set +u
if [[ -n "${study_name_base:-}" ]]; then
  study_name_base="${study_name_base//\'/}"
  base="$study_name_base"
elif [[ -n "${study_name:-}" ]]; then
  study_name="${study_name//\'/}"
  base="$study_name"
else
  echo "ERROR: Config must define study_name_base or study_name"
  exit 1
fi
set -u

sanitize() {
  echo "$1" | sed -E 's/[^A-Za-z0-9._-]+/_/g'
}

base_safe="$(sanitize "$base")"
model_safe="$(sanitize "$MODEL_KEY")"
metab_t_safe="$(sanitize "$METAB_TRANSFORM_KEY")"
taxa_t_safe="$(sanitize "$TAXA_TRANSFORM_KEY")"
metab_red_safe="$(sanitize "$metab_reduction")"
taxa_red_safe="$(sanitize "$taxa_reduction")"

study_name="${base_safe}__model-${model_safe}__metabT-${metab_t_safe}__taxaT-${taxa_t_safe}__metabR-${metab_red_safe}__taxaR-${taxa_red_safe}"
echo "Resolved study_name: $study_name"

module purge
module load miniforge3
unset PYTHONPATH PYTHONHOME

set +u
conda_sh="$(dirname "$(dirname "$(type -P conda)")")/etc/profile.d/conda.sh"
source "$conda_sh"
conda activate /work/PCDC/shared/IntegratedLearner_202412/conda_envs/R4
set -u

echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-unset}"
echo "SLURM_STEP_GPUS=${SLURM_STEP_GPUS:-unset}"
nvidia-smi -L || true
nvidia-smi --query-gpu=index,name,uuid --format=csv || true

echo "Using Rscript: $(which Rscript)"
Rscript --version || true
Rscript -e "library(xgboost); cat('xgboost=', as.character(packageVersion('xgboost')), '\n')"

cmd=(
  Rscript automated_pipeline.R
  --model_to_run "$model_to_run"
  --file_path "$file_path"
  --input_file "$input_file"
  --study_name "$study_name"
  --type_of_analysis "$type_of_analysis"
  --response_variable "$response_variable"
  --stratify_variable "$stratify_variable"
  --training_proportion "$training_proportion"
  --num_repeats "$num_repeats"
  --num_folds "$num_folds"
  --metab_transform "$metab_transform"
  --taxa_transform "$taxa_transform"
  --metab_reduction "$metab_reduction"
  --taxa_reduction "$taxa_reduction"
)

if [[ "$MODEL_KEY" == "xgb" ]]; then
  cmd+=(
    --use_gpu "$MODEL_USE_GPU"
    --gpu_device "$MODEL_GPU_DEVICE"
    --nthread "$LOCAL_CPUS_PER_TASK"
  )
fi

echo "Running command:"
printf ' %q' "${cmd[@]}"
echo

"${cmd[@]}"
