#!/bin/bash
#SBATCH --job-name=genome_assembly_and_QC
#SBATCH --output=Logs/genome_assembly_and_QC_%j.out
#SBATCH --error=Logs/genome_assembly_and_QC_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal

# Load conda
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

echo "-------------------------------------------------"
echo "Started genome assembly and QC job on SLURM"
echo "Using $SLURM_CPUS_PER_TASK CPUs"
echo "Job started at: $(date)"
echo "-------------------------------------------------"

# Run assembly pipeline
bash Pipeline.sh

echo "-------------------------------------------------"
echo "Job finished at: $(date)"
echo "-------------------------------------------------"
