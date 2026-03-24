#!/bin/bash
#SBATCH --job-name=dataset_anaysis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=128:00:00
#SBATCH --output=logs/dataset_analysis%j.out
#SBATCH --error=logs/dataset_analysis_%j.err

source ~/.bashrc
source ~/.nextflowrc
nextflow 00_dataset_analysis.nf -entry PROTEIN_RMSD_ANALYSIS -c ../nextflow.config -resume
