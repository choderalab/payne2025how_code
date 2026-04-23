#!/bin/bash
#SBATCH --job-name=plif_recall
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/plif_recall_%j.out
#SBATCH --error=logs/plif_recall_%j.err

source ~/.bashrc
source ~/.nextflowrc

nextflow 04_docking_analysis_v2.nf -entry PLIF_RECALL_POSIT -c ../nextflow.config -c nextflow.config -resume
