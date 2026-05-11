#!/bin/zsh
docked_sdf="/Users/apayne/Downloads/full_cross_dock_v2/docked_files/ALL_1_poses/142_docked/docking_results.sdf"
DATA="/Users/apayne/Downloads/full_cross_dock_v2"
REPO="/Users/apayne/science/sars-cov-2-retro-paper-asap"
SCRIPT="$REPO/nextflow_workflows/04_docking_analysis/scripts/calculate_plif_recall.py"
CACHE_DIR="$DATA/mpro_fragalysis-04-01-24_curated_cache_fixed"
CMPD_DICT="$DATA/cmpd_date_dict/structure_to_cmpd_dict.json"
python3 $SCRIPT \
        --docked-sdf "${docked_sdf}" \
        --cache-dir "${CACHE_DIR}" \
        --cmpd-dict "${CMPD_DICT}" \
        --output-csv "plif_recall_test.csv"