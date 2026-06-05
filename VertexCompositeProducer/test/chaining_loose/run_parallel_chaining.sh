#!/usr/bin/env bash
set -euo pipefail

MAX_I="${1:-23}"
MAX_JOBS="${MAX_JOBS:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs_parallel_chaining"

mkdir -p "${LOG_DIR}"

cd "${SCRIPT_DIR}"

echo "Running HIForward jobs 0..${MAX_I} with ${MAX_JOBS} parallel workers"

for i in $(seq 0 "${MAX_I}"); do
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        wait -n
    done

    echo "Starting i=${i}"
    root -l -b -q "sampleAna_chaining_by_i.C(${i})" > "${LOG_DIR}/chaining_${i}.log" 2>&1 &
done

wait
echo "All jobs finished. Logs are in ${LOG_DIR}"