#!/bin/bash
#SBATCH --partition=super
#SBATCH --mem=40G

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread

#SBATCH --mail-type=ALL
#SBATCH --mail-user=suzette.palmer@utsouthwestern.edu

#SBATCH --output=/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code/B_BashScripts/logs/%x_%j.out
#SBATCH --error=/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code/B_BashScripts/logs/%x_%j.err

set -euo pipefail

: "${CONFIG_SUBDIR:?ERROR: CONFIG_SUBDIR not set}"
: "${TASK_IDS:?ERROR: TASK_IDS not set}"

BASE_DIR="${BASE_DIR:-/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code}"
WORKER_SCRIPT="${BASE_DIR}/B_BashScripts/submit_integrative_array_small_single_node.sh"

if [[ ! -f "$WORKER_SCRIPT" ]]; then
  echo "ERROR: WORKER_SCRIPT not found: $WORKER_SCRIPT"
  exit 1
fi

read -r -a TASK_ID_LIST <<< "$TASK_IDS"

if (( ${#TASK_ID_LIST[@]} == 0 )); then
  echo "ERROR: No task IDs provided in TASK_IDS"
  exit 1
fi

LOCAL_NODE_NAME="${SLURMD_NODENAME:-$(hostname)}"

echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_NODELIST=${SLURM_JOB_NODELIST:-$LOCAL_NODE_NAME}"
echo "CONFIG_SUBDIR=${CONFIG_SUBDIR}"
echo "WORKER_SCRIPT=${WORKER_SCRIPT}"
echo "Running ${#TASK_ID_LIST[@]} task(s) on one node"
echo "Task IDs: ${TASK_ID_LIST[*]}"
echo

pids=()
task_ids_started=()

for task_id in "${TASK_ID_LIST[@]}"; do
  echo "[LAUNCH] task=${task_id}"

  srun \
    --exclusive \
    -N1 \
    -n1 \
    -c1 \
    --cpu-bind=cores \
    --export=ALL,CONFIG_SUBDIR="${CONFIG_SUBDIR}",BASE_DIR="${BASE_DIR}",SLURM_ARRAY_TASK_ID="${task_id}",LOCAL_CPUS_PER_TASK=1,LOCAL_JOB_ID="${SLURM_JOB_ID}_${task_id}",LOCAL_NODE_NAME="${LOCAL_NODE_NAME}",LOCAL_MEM="40G_shared" \
    bash "$WORKER_SCRIPT" &

  pids+=("$!")
  task_ids_started+=("$task_id")
done

rc=0
for i in "${!pids[@]}"; do
  pid="${pids[$i]}"
  tid="${task_ids_started[$i]}"
  if wait "$pid"; then
    echo "[DONE] task=${tid}"
  else
    echo "[FAIL] task=${tid}"
    rc=1
  fi
done

exit "$rc"
