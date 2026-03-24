#!/bin/bash
#SBATCH --job-name=collect_docking_results
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=128:00:00
#SBATCH --output=logs/collect_docking_results_%j.out
#SBATCH --error=logs/collect_docking_results_%j.err

source ~/.bashrc
source ~/.nextflowrc
nextflow 00_collect_docking_results.nf ../nextflow.config -resume
