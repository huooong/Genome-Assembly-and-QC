###################################
#!/bin/bash

# Input and output directory
INPUT_DIR="Data"
OUTPUT_DIR="Results"
FILTERED_DIR="$OUTPUT_DIR/Filtered"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Conda environments
FLYE="flye_env"
###################################

# Assembly with Flye
echo "---------------------------------------------------"
echo "Assembling sequenced data using Flye Assembler"
echo "---------------------------------------------------"

## Load Flye environment
conda activate $FLYE

## Directory for assemblies
ASSEMBLY_DIR="$OUTPUT_DIR/Assembly"
mkdir -p "$ASSEMBLY_DIR"

## Flye parameters
GENOME_SIZE="5m"
THREADS=4		#change

## Run flye
for file in "$FILTERED_DIR"/*/*_filtered.fastq; do
	sample=$(basename "$file" _filtered.fastq)
	echo "Assembling $sample..."

	flye \
		--nano-hq "$file" \
		--genome-size "$GENOME_SIZE" \
		--out-dir "$ASSEMBLY_DIR/$sample" \
		--threads "$THREADS"

echo "Assembly saved in $ASSEMBLY_DIR"
echo "---------------------------------------------------"

done
conda deactivate
