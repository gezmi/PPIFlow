#!/bin/bash
# Run AlphaFold 3 on a single JSON input.
#
# Usage:
#   bash workflow/run_af3_single.sh <json_path> <output_dir> <af3_home> [model_dir]
#
# Arguments:
#   json_path  - Path to AF3 JSON input file
#   output_dir - Directory for AF3 output
#   af3_home   - Path to AF3 installation (also used as conda env name)
#   model_dir  - (optional) Path to AF3 model weights directory

set -euo pipefail

JSON_PATH=$1
OUTPUT_DIR=$2
AF3_HOME=$3
MODEL_DIR=${4:-}

mkdir -p "$OUTPUT_DIR"

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "$AF3_HOME"

python3 "$AF3_HOME/alphafold3/run/run_af3.py" \
    --json_path "$JSON_PATH" \
    --output_dir "$OUTPUT_DIR" \
    ${MODEL_DIR:+--model_dir "$MODEL_DIR"} \
    --max_template_date 1900-01-01
