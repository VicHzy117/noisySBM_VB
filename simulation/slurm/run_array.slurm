#!/bin/bash -l
#SBATCH --job-name=paper-simulation
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=results/data/reruns/logs/%A_%a.out
#SBATCH --error=results/data/reruns/logs/%A_%a.err
set -euo pipefail
module load R/4.4.2-gfbf-2024a
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
Rscript simulation/run_simulation.R --task-id="${SLURM_ARRAY_TASK_ID}"
