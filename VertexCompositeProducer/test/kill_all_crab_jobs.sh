./#!/bin/bash

set -euo pipefail

shopt -s nullglob

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

mapfile -t crab_projects < <(
    find . -maxdepth 2 -type d -name 'crab_*' -path './crab_projects*/*' | sort
)

if [[ ${#crab_projects[@]} -eq 0 ]]; then
    echo "No CRAB project directories found under: $script_dir"
    exit 0
fi

echo "Found ${#crab_projects[@]} CRAB project directories. Killing all jobs..."

for project_dir in "${crab_projects[@]}"; do
    echo "Killing jobs in: $project_dir"
    crab kill -d "$project_dir"
done

echo "Done."