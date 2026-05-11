#!/bin/zsh
# Run a sample of PLIF recall jobs and extrapolate to the full 412-job run.
# Usage: ./time_plif_recall.sh [N_SAMPLE]   (default: 5)

N_SAMPLE=${1:-5}
TOTAL_JOBS=412

DATA="/Users/apayne/Downloads/full_cross_dock_v2"
REPO="/Users/apayne/science/sars-cov-2-retro-paper-asap"
SCRIPT="$REPO/nextflow_workflows/04_docking_analysis/scripts/calculate_plif_recall.py"
CACHE_DIR="$DATA/mpro_fragalysis-04-01-24_curated_cache_fixed"
CMPD_DICT="$DATA/cmpd_date_dict/structure_to_cmpd_dict.json"
SDF_DIR="$DATA/docked_files/ALL_1_poses"

# Pick N_SAMPLE jobs spread across the full list
sdfs=("$SDF_DIR"/*/docking_results.sdf)
n_total=${#sdfs[@]}
step=$(( n_total / N_SAMPLE ))

echo "Timing $N_SAMPLE jobs (every ~${step}th of $n_total)..."
echo "------------------------------------------------------------"

total_seconds=0
i=0
job_num=0
for sdf in "${sdfs[@]}"; do
    if (( i % step == 0 )) && (( job_num < N_SAMPLE )); then
        job_num=$(( job_num + 1 ))
        job_name=$(basename $(dirname $sdf))
        echo -n "Job $job_num/$N_SAMPLE ($job_name)... "
        t0=$SECONDS
        conda run -n drugforge python3 "$SCRIPT" \
            --docked-sdf "$sdf" \
            --cache-dir "$CACHE_DIR" \
            --cmpd-dict "$CMPD_DICT" \
            --output-csv "/tmp/plif_timing_${job_name}.csv" 2>/dev/null
        elapsed=$(( SECONDS - t0 ))
        total_seconds=$(( total_seconds + elapsed ))
        echo "${elapsed}s"
    fi
    i=$(( i + 1 ))
done

avg=$(( total_seconds / N_SAMPLE ))
total_est=$(( avg * TOTAL_JOBS ))
total_min=$(( total_est / 60 ))
total_hr=$(echo "scale=1; $total_est / 3600" | bc)

echo "------------------------------------------------------------"
echo "Mean per job:      ${avg}s"
echo "Estimated total:   ${total_min} min (${total_hr} hr) for $TOTAL_JOBS jobs"
echo "  (assumes serial; on cluster these run in parallel)"
