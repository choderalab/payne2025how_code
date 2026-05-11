#!/bin/bash
#SBATCH --job-name=plif_top4_filtered
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/plif_top4_filtered_%j.out
#SBATCH --error=logs/plif_top4_filtered_%j.err

source ~/.bashrc
source ~/.nextflowrc

# Filters ALL_1_poses_plif.parquet to top-4 scaffolds only, then evaluates:
#   self_docked:       top-4 ref == top-4 query scaffold
#   not_top4_to_top4:  top-4 ref, non-top-4 query scaffold
# Both with PLIF Tversky recall cutoffs 0.5/0.75/1.0 and RMSD+POSIT scorers.
nextflow 04_docking_analysis_v2.nf -entry analyze_posit_plif_top4_filtered -c ../nextflow.config -c nextflow.config -resume
