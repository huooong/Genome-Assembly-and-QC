###################################
#!/bin/bash

# Input and output directory
INPUT_DIR="Data"
OUTPUT_DIR="Results"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

###################################

# Filter reads using Filtlong
echo "---------------------------------------------------"
echo "Filtering reads using Filtlong"
echo "---------------------------------------------------"

## Load Filtlong environment
conda activate $FILTLONG

## Directory for filtered results
FILTERED_DIR="$OUTPUT_DIR/Filtered"
mkdir -p "$FILTERED_DIR"

## Filtlong parameters 
MIN_LENGTH=1000        #Minimum read length to keep
KEEP_PERCENT=90        #Percentage of best reads to retain
TARGET_BASES=500000000   

## Run Filtlong
for file in "$INPUT_DIR"/*.fastq.gz; do
    sample=$(basename "$file" .fastq.gz)
    sample_dir="$FILTERED_DIR/$sample"
    mkdir -p "$sample_dir"
    
    echo "Filtering $file ....."

    filtlong \
        --min_length "$MIN_LENGTH" \
        --keep_percent "$KEEP_PERCENT" \
        --target_bases "$TARGET_BASES" \
        "$file" > "$sample_dir/${sample}_filtered.fastq"

    echo "Filtered $sample has been saved in $sample_dir"
    echo "---------------------------------------------------"

done
conda deactivate
