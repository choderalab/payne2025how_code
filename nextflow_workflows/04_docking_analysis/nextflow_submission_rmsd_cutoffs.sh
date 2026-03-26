#!/bin/bash
#SBATCH --job-name=rmsd_cutoff_sensitivity
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=128:00:00
#SBATCH --output=logs/rmsd_cutoff_sensitivity_%j.out
#SBATCH --error=logs/rmsd_cutoff_sensitivity_%j.err

source ~/.bashrc
source ~/.nextflowrc

# Step 1: generate evaluator factory settings for 1.5 and 2.5 Å cutoffs
nextflow 04_docking_analysis_v2.nf -entry GENERATE_SETTINGS_RMSD_1P5 -c nextflow.config -resume
nextflow 04_docking_analysis_v2.nf -entry GENERATE_SETTINGS_RMSD_2P5 -c nextflow.config -resume

# Step 2: run the full POSIT analysis at each cutoff
nextflow 04_docking_analysis_v2.nf -entry analyze_posit_rmsd_1p5 -c nextflow.config -resume
nextflow 04_docking_analysis_v2.nf -entry analyze_posit_rmsd_2p5 -c nextflow.config -resume
