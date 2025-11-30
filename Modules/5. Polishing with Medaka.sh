###################################
#!/bin/bash

# Directories
INPUT_DIR="Data"
OUTPUT_DIR="Results"
FILTERED_DIR="$OUTPUT_DIR/Filtered"
ASSEMBLY_DIR="$OUTPUT_DIR/Assembly"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Parameters
RACON_ROUNDS=3 # Change if needed

###################################

# Final polishing with Medaka
echo "---------------------------------------------------"
echo "Final polishing with Medaka"
echo "---------------------------------------------------"

## Medaka parameters
THREADS=8        # I changed this to 16 when I ran Medaka
MODEL="r1041_e82_400bps_hac_v5.0.0"        

## Directory for Medaka
MEDAKA_DIR="$OUTPUT_DIR/Medaka"
mkdir -p "$MEDAKA_DIR"

## Load Medaka environment
conda activate medaka

## Run Medaka on Racon-polished assemblies
for polished_assembly in "$RACON_DIR"/*; do
    sample=$(basename "$polished_assembly")
    sample_dir="$MEDAKA_DIR/$sample"
    mkdir -p "$sample_dir"

    echo "Starting polishing with Medaka for $sample ....."

    racon_polished_assembly="$RACON_DIR/${sample}_racon_round_${RACON_ROUNDS}.fasta"
    reads="$FILTERED_DIR/$sample/${sample}_filtered.fastq"
    
    medaka_consensus \
        -i $reads \
        -d $racon_polished_assembly \
        -o $sample_dir \
        -t $THREADS \
        -m $MODEL

    echo
    echo "Completed Medaka for $sample. Results have been saved in "$MEDAKA_DIR""
    echo "---------------------------------------------------"
done
conda deactivate

