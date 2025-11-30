###################################
#!/bin/bash

# Directories
INPUT_DIR="Data"
OUTPUT_DIR="Results"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Conda environment
NANOPLOT="nanoplot2"
###################################

# Assess read quality of raw reads with NanoPlot
echo "---------------------------------------------------"
echo "Assessing read quality of raw reads with NanoPlot"
echo "---------------------------------------------------"

## Load NanoPlot environment
conda activate $NANOPLOT

## Directory for NanoPlot Results
NANOPLOT_DIR="$OUTPUT_DIR/NanoPlot"
mkdir -p "$NANOPLOT_DIR"

## Run NanoPlot
for file in "$INPUT_DIR"/*.fastq.gz; do
    sample=$(basename "$file" .fastq.gz)
    sample_dir="$NANOPLOT_DIR/$sample"
    mkdir -p "$sample_dir"

    echo "Running NanoPlot on $sample ....."
    NanoPlot \
        --fastq "$file" \
        -o "$sample_dir" \
        --prefix "$sample" \
        --tsv_stats --raw
    
    echo 
    
    echo "NanoPlot results for $sample have been saved in $sample_dir"
    echo "---------------------------------------------------"

done
conda deactivate
