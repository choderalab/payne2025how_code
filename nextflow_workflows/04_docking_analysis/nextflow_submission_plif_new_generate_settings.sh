#!/bin/bash
#SBATCH --job-name=plif_generate_settings
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/plif_generate_settings_%j.out
#SBATCH --error=logs/plif_generate_settings_%j.err

source ~/.bashrc
source ~/.nextflowrc

# PLIF Tversky recall cutoffs 0.75 and 1.0 (generates all split YAMLs for each)
nextflow 04_docking_analysis_v2.nf -entry GENERATE_SETTINGS_PLIF_0_75 -c ../nextflow.config -c nextflow.config
nextflow 04_docking_analysis_v2.nf -entry GENERATE_SETTINGS_PLIF_1_0  -c ../nextflow.config -c nextflow.config

# PLIF Recall as scoring method (generates 4 YAMLs: rmsd2, plif0.5, plif0.75, plif1.0)
nextflow 04_docking_analysis_v2.nf -entry GENERATE_SETTINGS_PLIF_SCORER -c ../nextflow.config -c nextflow.config
