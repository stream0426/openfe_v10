#!/bin/bash
#SBATCH --job-name=openfe-rbfe
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

set -euo pipefail

SIF_PATH=${SIF_PATH:-/path/to/openfe.sif}
PROJECT_DIR=${PROJECT_DIR:-$PWD}
CONFIG_PATH=${CONFIG_PATH:-$PROJECT_DIR/configs/rbfe_config.example.yaml}
TRANSFORM_JSON=${TRANSFORM_JSON:-$PROJECT_DIR/outputs/tyk2_rbfe/transformations/example.json}
SCRATCH_DIR=${SCRATCH_DIR:-${SLURM_TMPDIR:-/tmp}/openfe-${SLURM_JOB_ID:-manual}}

mkdir -p "$SCRATCH_DIR" "$PROJECT_DIR/logs"

apptainer exec --nv \
  --bind "$PROJECT_DIR":"$PROJECT_DIR" \
  --bind "$SCRATCH_DIR":/tmp \
  "$SIF_PATH" \
  /opt/conda/envs/openfe/bin/python \
  "$PROJECT_DIR/scripts/build_rbfe_network.py" \
  "$CONFIG_PATH"

apptainer exec --nv \
  --bind "$PROJECT_DIR":"$PROJECT_DIR" \
  --bind "$SCRATCH_DIR":/tmp \
  "$SIF_PATH" \
  /opt/conda/envs/openfe/bin/openfe \
  quickrun \
  "$TRANSFORM_JSON" \
  -d "$PROJECT_DIR/run-$(basename "${TRANSFORM_JSON%.json}")"

