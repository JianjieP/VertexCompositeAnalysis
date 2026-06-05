#!/usr/bin/env bash

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_dir="$script_dir"
project_glob="crab_projects*"
task_glob="*"
log_name="crab_status.log"
dry_run=false
use_existing_log=false
summary_file=""
max_resubmits=0

status_names=(
    finished
    failed
    idle
    running
    transferring
    unsubmitted
    cooloff
    held
    killed
    unknown
)

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Check CRAB task status under crab_projects*/* and resubmit failed jobs.

Options:
  -b, --base-dir DIR        Directory containing crab project folders.
                            Default: directory of this script
  -p, --project-glob GLOB   Project directory glob under base-dir.
                            Default: crab_projects*
  -t, --task-glob GLOB      Task directory glob under project dirs.
                            Default: *
  -n, --dry-run             Check and report only; do not run crab resubmit.
  -e, --use-existing-log    Do not run crab status; parse existing logs.
                            Looks for LOG_NAME first, then crab.log.
  -l, --log-name NAME       Status stdout log name in each task directory.
                            Default: crab_status.log
  -s, --summary-file FILE   TSV summary output path.
                            Default: BASE_DIR/check_and_resubmit_summary_YYYYmmdd_HHMMSS.tsv
  -m, --max-resubmits N     Stop resubmitting after N tasks. 0 means no limit.
  -h, --help                Show this help message.

Examples:
  $(basename "$0") --dry-run
  $(basename "$0") --project-glob 'crab_projects_0520' --task-glob 'crab_HIForward[0-9]_*'
  $(basename "$0") --use-existing-log --dry-run
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

warn() {
    echo "WARNING: $*" >&2
}

is_nonnegative_int() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--base-dir)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            base_dir="$2"
            shift 2
            ;;
        -p|--project-glob)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            project_glob="$2"
            shift 2
            ;;
        -t|--task-glob)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            task_glob="$2"
            shift 2
            ;;
        -n|--dry-run)
            dry_run=true
            shift
            ;;
        -e|--use-existing-log)
            use_existing_log=true
            shift
            ;;
        -l|--log-name)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            log_name="$2"
            shift 2
            ;;
        -s|--summary-file)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            summary_file="$2"
            shift 2
            ;;
        -m|--max-resubmits)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            is_nonnegative_int "$2" || die "--max-resubmits must be a non-negative integer"
            max_resubmits="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -d "$base_dir" ]] || die "Base directory does not exist: $base_dir"

if [[ "$use_existing_log" == false ]]; then
    command -v crab >/dev/null 2>&1 || die "Cannot find 'crab' in PATH. Set up the CRAB environment first, or use --use-existing-log."
fi

if [[ -z "$summary_file" ]]; then
    summary_file="$base_dir/check_and_resubmit_summary_$(date +%Y%m%d_%H%M%S).tsv"
fi

extract_status_count() {
    local state="$1"
    local log_file="$2"

    sed -nE "s/^(Jobs status:[[:space:]]*)?[[:space:]]*${state}[[:space:]]+[<>]?[0-9.]+%[[:space:]]*\\([[:space:]]*([0-9]+)\\/([0-9]+)\\).*/\\2/p" "$log_file" | head -n 1
}

extract_total_jobs() {
    local log_file="$1"

    sed -nE "s/^(Jobs status:[[:space:]]*)?[[:space:]]*(failed|finished|idle|running|transferring|unsubmitted|cooloff|held|killed|unknown)[[:space:]]+[<>]?[0-9.]+%[[:space:]]*\\([[:space:]]*[0-9]+\\/([0-9]+)\\).*/\\3/p" "$log_file" | head -n 1
}

status_count_or_zero() {
    local state="$1"
    local log_file="$2"
    local count

    count="$(extract_status_count "$state" "$log_file")"
    if [[ "$count" =~ ^[0-9]+$ ]]; then
        echo "$count"
    else
        echo 0
    fi
}

pick_existing_log() {
    local task_dir="$1"
    local preferred_log="$task_dir/$log_name"
    local crab_log="$task_dir/crab.log"

    if [[ -s "$preferred_log" ]]; then
        echo "$preferred_log"
    elif [[ -s "$crab_log" ]]; then
        echo "$crab_log"
    else
        return 1
    fi
}

write_summary_header() {
    local file="$1"
    printf 'task_dir\tstatus_ok\tfinished\tfailed\tidle\trunning\ttransferring\tunsubmitted\tcooloff\theld\tkilled\tunknown\ttotal\tunfinished\tresubmit_action\n' > "$file"
}

append_summary_row() {
    local file="$1"
    local task_dir="$2"
    local status_ok="$3"
    local total="$4"
    local unfinished="$5"
    local action="$6"
    shift 6
    local counts=("$@")

    printf '%s\t%s' "$task_dir" "$status_ok" >> "$file"
    for count in "${counts[@]}"; do
        printf '\t%s' "$count" >> "$file"
    done
    printf '\t%s\t%s\t%s\n' "$total" "$unfinished" "$action" >> "$file"
}

print_banner() {
    local message="$1"
    echo "##############################################"
    echo "## $message"
    echo "##############################################"
}

task_pattern="$base_dir/$project_glob/$task_glob"
mapfile -t task_dirs < <(compgen -G "$task_pattern")

filtered_task_dirs=()
for task_dir in "${task_dirs[@]}"; do
    [[ -d "$task_dir" ]] && filtered_task_dirs+=( "$task_dir" )
done
task_dirs=( "${filtered_task_dirs[@]}" )

(( ${#task_dirs[@]} > 0 )) || die "No CRAB task directories matched: $task_pattern"

write_summary_header "$summary_file" || die "Cannot write summary file: $summary_file"

total_tasks=0
status_ok_tasks=0
status_failed_tasks=0
parse_failed_tasks=0
finished_tasks=0
unfinished_jobs_total=0
total_jobs_all=0
resubmission_count=0
resubmit_failed_count=0
dry_run_resubmit_count=0

declare -A totals_by_status=()
for status in "${status_names[@]}"; do
    totals_by_status["$status"]=0
done

failed_tasks=()
unfinished_tasks=()
status_problem_tasks=()

for task_dir in "${task_dirs[@]}"; do
    ((total_tasks++))

    print_banner "Checking task: $task_dir"

    status_log="$task_dir/$log_name"
    status_ok=true
    resubmit_action="none"

    if [[ "$use_existing_log" == true ]]; then
        if ! status_log="$(pick_existing_log "$task_dir")"; then
            warn "No existing $log_name or crab.log found for $task_dir"
            ((status_failed_tasks++))
            status_problem_tasks+=( "$task_dir" )
            append_summary_row "$summary_file" "$task_dir" "no_log" 0 0 "none" 0 0 0 0 0 0 0 0 0 0
            echo
            continue
        fi
        echo "Using existing log: $status_log"
        ((status_ok_tasks++))
    else
        echo "Writing status output to: $status_log"
        if ! crab status -d "$task_dir" 2>&1 | tee "$status_log"; then
            warn "crab status failed for $task_dir"
            status_ok=false
            ((status_failed_tasks++))
            status_problem_tasks+=( "$task_dir" )
        else
            ((status_ok_tasks++))
        fi
    fi

    counts=()
    for status in "${status_names[@]}"; do
        count="$(status_count_or_zero "$status" "$status_log")"
        counts+=( "$count" )
        ((totals_by_status["$status"]+=count))
    done

    finished_count="${counts[0]}"
    failed_count="${counts[1]}"
    total_jobs="$(extract_total_jobs "$status_log")"

    if ! [[ "$total_jobs" =~ ^[0-9]+$ ]]; then
        warn "Could not parse total job count for $task_dir"
        ((parse_failed_tasks++))
        status_problem_tasks+=( "$task_dir" )
        total_jobs=0
        unfinished_jobs=0
    else
        unfinished_jobs=$((total_jobs - finished_count))
        ((unfinished_jobs < 0)) && unfinished_jobs=0
        ((total_jobs_all+=total_jobs))
        ((unfinished_jobs_total+=unfinished_jobs))

        if (( total_jobs > 0 && finished_count == total_jobs )); then
            ((finished_tasks++))
        elif (( unfinished_jobs > 0 )); then
            unfinished_tasks+=( "$task_dir ($finished_count/$total_jobs finished)" )
        fi
    fi

    if (( failed_count > 0 )); then
        failed_tasks+=( "$task_dir ($failed_count failed)" )

        if [[ "$dry_run" == true ]]; then
            echo ">>> DRY-RUN: would resubmit $failed_count failed job(s) in $task_dir"
            resubmit_action="dry-run"
            ((dry_run_resubmit_count++))
        elif (( max_resubmits > 0 && resubmission_count >= max_resubmits )); then
            echo ">>> Skipping resubmit because --max-resubmits=$max_resubmits has been reached."
            resubmit_action="max-skip"
        else
            echo ">>> Resubmitting $failed_count failed job(s) in $task_dir"
            if crab resubmit -d "$task_dir"; then
                resubmit_action="resubmitted"
                ((resubmission_count++))
            else
                warn "crab resubmit failed for $task_dir"
                resubmit_action="resubmit-failed"
                ((resubmit_failed_count++))
            fi
        fi
    fi

    append_summary_row "$summary_file" "$task_dir" "$status_ok" "$total_jobs" "$unfinished_jobs" "$resubmit_action" "${counts[@]}"

    echo "Status: finished=$finished_count/$total_jobs failed=$failed_count unfinished=$unfinished_jobs"
    echo "## Done with task: $task_dir"
    echo
done

print_list() {
    local title="$1"
    shift
    local items=("$@")

    [[ ${#items[@]} -gt 0 ]] || return 0
    echo
    echo "$title"
    printf '  %s\n' "${items[@]}"
}

print_banner "Summary"
echo "Base directory:           $base_dir"
echo "Task pattern:             $project_glob/$task_glob"
echo "Total tasks found:        $total_tasks"
echo "Status OK tasks:          $status_ok_tasks"
echo "Status/log problem tasks: $status_failed_tasks"
echo "Parse problem tasks:      $parse_failed_tasks"
echo "Fully finished tasks:     $finished_tasks"
echo "Unfinished jobs:          $unfinished_jobs_total / $total_jobs_all"
echo "Resubmissions made:       $resubmission_count"
echo "Resubmissions failed:     $resubmit_failed_count"
echo "Dry-run resubmits:        $dry_run_resubmit_count"
echo "Summary TSV:              $summary_file"
echo
echo "OVERALL JOB STATUS ACROSS ALL TASKS:"
for status in "${status_names[@]}"; do
    printf '%-14s %s\n' "$status:" "${totals_by_status[$status]}"
done
echo "TOTAL JOBS:     $total_jobs_all"

print_list "Tasks with failed jobs:" "${failed_tasks[@]}"
print_list "Tasks not fully finished:" "${unfinished_tasks[@]}"
print_list "Tasks with status/log/parse problems:" "${status_problem_tasks[@]}"

echo "##############################################"
