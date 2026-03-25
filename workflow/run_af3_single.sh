#!/bin/bash
# Run AlphaFold 3 on a single JSON input (no MSA generation).
#
# Usage:
#   bash workflow/run_af3_single.sh <json_path> <output_dir> <af3_home> [model_dir]
#
# Arguments:
#   json_path  - Path to AF3 JSON input file (must have empty pairedMsa/unpairedMsa)
#   output_dir - Directory for AF3 output
#   af3_home   - Path to AF3 installation (also used as conda env name)
#   model_dir  - (optional) Path to AF3 model weights directory

set -euo pipefail

JSON_PATH=$1
OUTPUT_DIR=$2
AF3_HOME=$3
MODEL_DIR=${4:-/home/labs/schreiber/$USER/models}

mkdir -p "$OUTPUT_DIR"

set +u  # conda activate scripts may reference unbound vars (e.g. ADDR2LINE from binutils)
eval "$('/apps/easybd/programs/miniconda/24.11_environmentally/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
conda activate "$AF3_HOME"
set -u
export PYTHONNOUSERSITE=1
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false --xla_gpu_graph_level=0"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.90

echo "Working on $JSON_PATH"
time python3 "$AF3_HOME/alphafold3/run/run_af3.py" \
    --json_path "$JSON_PATH" \
    --output_dir "$OUTPUT_DIR" \
    --model_dir "$MODEL_DIR" \
    --max_template_date 1900-01-01 \
    --flash_attention_implementation xla

echo "Job finished for $JSON_PATH"
