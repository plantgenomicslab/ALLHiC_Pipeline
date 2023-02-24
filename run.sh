#!/usr/bin/env bash

set -x
set -e

snakemake --rerun-incomplete \
          --printshellcmds \
          --max-jobs-per-second 5 \
          --max-status-checks-per-second 3 \
          --jobs 48 \
          --latency-wait 40 \
          --notemp \
          --cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} --parsable -A {cluster.account} -p {cluster.partition} -t {cluster.time} -o {cluster.output} -e {cluster.error}" \
          --cluster-config cluster.json \
          --cluster-status scripts/status.py \
          --cluster-cancel 'scancel' \
          "$@"


