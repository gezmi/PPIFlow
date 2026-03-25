#!/bin/bash
# LSF submission wrapper for Snakemake 7
# Routes GPU jobs to short-gpu queue, CPU jobs to short queue
#
# Called by Snakemake as:
#   lsf_submit.sh {threads} {resources.mem_mb} {resources.gpu} {rule} /path/to/jobscript.sh

THREADS=$1
MEM_MB=$2
GPU=$3
RULE=$4
shift 4
# Remaining argument(s): the job script path

QUEUE="short"
GPU_FLAG=""
if [ "$GPU" -ge 1 ] 2>/dev/null; then
    QUEUE="short-gpu"
    GPU_FLAG='-gpu "num=1"'
fi

LOGDIR="${SM_LOGDIR:-logs}"
mkdir -p "$LOGDIR"

exec bsub -q "$QUEUE" \
     -n "$THREADS" \
     -R "rusage[mem=${MEM_MB}]" \
     $GPU_FLAG \
     -J "sm_${RULE}" \
     -o "${LOGDIR}/sm_${RULE}_%J.out" \
     -e "${LOGDIR}/sm_${RULE}_%J.err" \
     "$@"
