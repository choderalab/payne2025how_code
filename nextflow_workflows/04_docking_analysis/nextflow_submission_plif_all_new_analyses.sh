#!/bin/bash
#SBATCH --job-name=plif_all_new_analyses
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/plif_all_new_analyses_%j.out
#SBATCH --error=logs/plif_all_new_analyses_%j.err

source ~/.bashrc
source ~/.nextflowrc

# Runs all three new PLIF analyses concurrently within a single Nextflow invocation:
#   analyze_posit_plif_cutoffs      — PLIF recall at 0.5/0.75/1.0 with datesplit + scaffold splits
#   analyze_posit_plif_scorer       — PLIF recall as pose-selection scorer, 4 eval metrics
#   analyze_posit_plif_top4_filtered — top-4 scaffold self-docked and not-top4-to-top4
#
# Run nextflow_submission_plif_new_generate_settings.sh first to generate the required YAMLs.
nextflow 04_docking_analysis_v2.nf -entry analyze_all_new_plif -c ../nextflow.config -c nextflow.config -resume
