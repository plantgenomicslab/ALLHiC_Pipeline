#!/usr/bin/env bash

set -x
set -e

snakemake --rerun-incomplete --cluster-config cluster.json \
          --cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} --parsable -A {cluster.account} -p {cluster.partition} -t {cluster.time} -o {cluster.output} -e {cluster.error}" \
          --max-jobs-per-second 5 \
          --max-status-checks-per-second 3 \
		  --jobs 4 \
		  --latency-wait 40 \
		  --notemp \
          "$@"
