#!/bin/bash
# Snakemake LSF cluster submission wrapper.
#
# Usage:
#   snakemake --cluster "bash workflow/lsf_submit.sh {resources.gpu} \
#     {resources.mem_mb} {threads} {rule}" ...
#
# Handles GPU vs CPU-only jobs automatically.

GPU=$1
MEM_MB=$2
THREADS=$3
RULE=$4
shift 4

mkdir -p logs

GPU_FLAG=""
QUEUE="short"
if [ "$GPU" -gt 0 ] 2>/dev/null; then
    GPU_FLAG="-gpu num=${GPU}"
    QUEUE="short-gpu"
fi

bsub -q "$QUEUE" \
    $GPU_FLAG \
    -n "$THREADS" \
    -R "rusage[mem=${MEM_MB}]" \
    -W 6:00 \
    -J "$RULE" \
    -o "logs/${RULE}.%J.out" \
    -e "logs/${RULE}.%J.err" \
    "$@"
