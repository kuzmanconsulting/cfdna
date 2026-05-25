#!/bin/bash
set -e

nohup micromamba run -n snakemake \
    snakemake --profile profiles/default > snakemake.log 2>&1 & disown
echo $! > snakemake.pid
echo "Launched snakemake (PID $(cat snakemake.pid)) — logs: snakemake.log"
