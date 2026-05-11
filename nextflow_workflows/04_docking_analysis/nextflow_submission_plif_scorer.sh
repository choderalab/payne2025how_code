#!/bin/bash
#SBATCH --job-name=plif_scorer_analysis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/plif_scorer_analysis_%j.out
#SBATCH --error=logs/plif_scorer_analysis_%j.err

source ~/.bashrc
source ~/.nextflowrc

# Select poses by highest PLIF recall; evaluate against RMSD<2, PLIF>=0.5/0.75/1.0
nextflow 04_docking_analysis_v2.nf -entry analyze_posit_plif_scorer -c ../nextflow.config -c nextflow.config -resume
