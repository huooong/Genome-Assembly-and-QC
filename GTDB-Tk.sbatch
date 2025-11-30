#!/bin/bash
#SBATCH --job-name=GTDB-Tk
#SBATCH --output=logs/GTDB-Tk_%j.out
#SBATCH --error=logs/GTDB-Tk_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=150G
#SBATCH --partition=normal

# Load conda
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

echo "-------------------------------------------------"
echo "Started taxonomy asssignment job on SLURM"
echo "Using $SLURM_CPUS_PER_TASK CPUs"
echo "Job started at: $(date)"
echo "-------------------------------------------------"

# Run assembly pipeline
bash GTDB-Tk.sh

echo "-------------------------------------------------"
echo "Job finished at: $(date)"
echo "-------------------------------------------------"
