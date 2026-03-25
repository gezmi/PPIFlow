#!/bin/bash
# Snakemake cluster status script for LSF.
# Usage: --cluster-status "bash workflow/cluster_status.sh"
#
# Snakemake calls this with the job ID (from bsub output) to check status.
# Must print one of: running, success, failed

# Extract numeric job ID from bsub output like "Job <12345> is submitted..."
JOBID=$(echo "$1" | grep -oP '\d+' | head -1)

if [ -z "$JOBID" ]; then
    echo "failed"
    exit 0
fi

# Query job status
STATUS=$(bjobs -noheader -o "stat" "$JOBID" 2>/dev/null)

case "$STATUS" in
    PEND|RUN|PSUSP|USUSP|SSUSP|WAIT)
        echo "running"
        ;;
    DONE)
        echo "success"
        ;;
    EXIT|ZOMBI|"")
        echo "failed"
        ;;
    *)
        echo "failed"
        ;;
esac
