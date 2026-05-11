#!/bin/bash
#SBATCH --job-name=plif_cutoff_analysis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/plif_cutoff_analysis_%j.out
#SBATCH --error=logs/plif_cutoff_analysis_%j.err

source ~/.bashrc
source ~/.nextflowrc

# Datesplit + scaffold splits (x_to_x, not_x_to_x, x_to_not_x) at all three
# PLIF Tversky recall cutoffs. The 0.5 results are cached via -resume.
nextflow 04_docking_analysis_v2.nf -entry analyze_posit_plif_cutoffs -c ../nextflow.config -c nextflow.config -resume
