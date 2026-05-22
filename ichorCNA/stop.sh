#!/bin/bash
set -e

if [[ ! -f snakemake.pid ]]; then
    echo "No snakemake.pid found"
    exit 1
fi

PID=$(cat snakemake.pid)
if kill -0 "$PID" 2>/dev/null; then
    kill "$PID"
    echo "Sent SIGTERM to snakemake (PID $PID)"
else
    echo "PID $PID not running"
fi
