#./check_and_resubmit.sh

# 进入工作目录（换成你自己的）
#cd crab_projects

# Initialize counters for resubmissions, finished tasks, unfinished jobs, and total jobs
resubmission_count=0
finished_tasks_count=0
unfinished_jobs_count=0
total_jobs_count=0
total_count=0

# Initialize counters for total job status across all tasks
total_finished_jobs=0
total_idle_jobs=0
total_running_jobs=0
total_transferring_jobs=0
total_unsubmitted_jobs=0
total_failed_jobs=0

# Extract numerator (x in x/y) from a Jobs status line for a given state.
extract_status_count() {
    local state="$1"
    local log_file="$2"
    local line
    line=$(grep -E "^[[:space:]]+${state}[[:space:]]+[0-9.]+%[[:space:]]*\([[:space:]]*[0-9]+/[0-9]+\)" "$log_file" | head -1)
    if [[ -n "$line" ]]; then
        echo "$line" | sed -E 's/.*\(([[:space:]]*[0-9]+)\/[0-9]+\).*/\1/' | tr -d ' '
    else
        echo 0
    fi
}

# Extract denominator (y in x/y) from the first Jobs status line found.
extract_total_jobs() {
    local log_file="$1"
    local line
    line=$(grep -E "^[[:space:]]+(failed|finished|running|idle|transferring|unsubmitted)[[:space:]]+[0-9.]+%[[:space:]]*\([[:space:]]*[0-9]+/[0-9]+\)" "$log_file" | head -1)
    if [[ -n "$line" ]]; then
        echo "$line" | sed -E 's/.*\([[:space:]]*[0-9]+\/([0-9]+)\).*/\1/' | tr -d ' '
    else
        echo 0
    fi
}

# 遍历所有任务目录
for d in crab_projects/*; do
    echo "##############################################"
    echo "## Checking task: $d"
    echo "##############################################"
    
    # 保存状态输出
    crab status -d "$d" | tee "$d/crab.log"

    # Increment the total count
    ((total_count++))

    # If the task has any failed jobs
    if grep -q "Jobs status:                    failed" "$d/crab.log"; then
        echo ">>> Resubmitting failed jobs in $d"
        crab resubmit -d "$d"
        
        # Increment the resubmission counter
        ((resubmission_count++))
    fi

    # Extract finished jobs and total jobs from strict "Jobs status" lines only.
    finished_jobs=$(extract_status_count "finished" "$d/crab.log")
    total_jobs=$(extract_total_jobs "$d/crab.log")

    # Debugging: Print the extracted values
    echo "Debug: finished_jobs = $finished_jobs"
    echo "Debug: total_jobs = $total_jobs"

    # Extract individual job status counts for overall statistics
    finished_count=$(extract_status_count "finished" "$d/crab.log")
    idle_count=$(extract_status_count "idle" "$d/crab.log")
    running_count=$(extract_status_count "running" "$d/crab.log")
    transferring_count=$(extract_status_count "transferring" "$d/crab.log")
    unsubmitted_count=$(extract_status_count "unsubmitted" "$d/crab.log")
    failed_count=$(extract_status_count "failed" "$d/crab.log")

    # Add to totals (use 0 if empty)
    ((total_finished_jobs+=${finished_count:-0}))
    ((total_idle_jobs+=${idle_count:-0}))
    ((total_running_jobs+=${running_count:-0}))
    ((total_transferring_jobs+=${transferring_count:-0}))
    ((total_unsubmitted_jobs+=${unsubmitted_count:-0}))
    ((total_failed_jobs+=${failed_count:-0}))

    # Check if finished_jobs and total_jobs are valid numbers
    if [[ -n "$finished_jobs" && "$finished_jobs" =~ ^[0-9]+$ && -n "$total_jobs" && "$total_jobs" =~ ^[0-9]+$ ]]; then
        # Fully finished task
        if [ "$finished_jobs" -eq "$total_jobs" ]; then
            ((finished_tasks_count++))
        fi

        # Calculate unfinished jobs
        unfinished_jobs=$((total_jobs - finished_jobs))
        ((unfinished_jobs_count+=unfinished_jobs))
        ((total_jobs_count+=total_jobs))
    else
        echo "Error: Invalid data in task $d, skipping job status."
    fi

    echo "## Done with task: $d"
    echo
done

# Report the total number of tasks, fully finished tasks, unfinished jobs, and resubmissions
echo "##############################################"
echo "Total tasks processed: $total_count"
echo "Fully finished tasks (100%): $finished_tasks_count"
echo "Unfinished jobs: $unfinished_jobs_count / $total_jobs_count"
echo "Resubmissions made: $resubmission_count"
echo ""
echo "OVERALL JOB STATUS ACROSS ALL TASKS:"
echo "finished:      $total_finished_jobs"
echo "idle:          $total_idle_jobs"
echo "running:       $total_running_jobs"
echo "transferring:  $total_transferring_jobs"
echo "unsubmitted:   $total_unsubmitted_jobs"
echo "failed:        $total_failed_jobs"
echo "TOTAL JOBS:    $total_jobs_count"
echo "##############################################"
