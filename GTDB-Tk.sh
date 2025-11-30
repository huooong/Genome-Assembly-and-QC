###################################
#!/bin/bash

# Directories
INPUT_DIR="Data"
OUTPUT_DIR="Results"
MEDAKA_DIR="$OUTPUT_DIR/Medaka"
GENOME_DIR="$OUTPUT_DIR/CheckM/Genomes"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"
GTDBTK="gtdbtk-2.5.2"

###################################

# Assigning taxonomy with GTDB-Tk
echo "---------------------------------------------------"
echo "Assigning taxonomy with GTDB-Tk"
echo "---------------------------------------------------"

## Load GTDB-Tk environment
conda activate $GTDBTK

## GTDB-Tk parameters
GTDBTK_DATA_PATH="/databases/GTDB/gtdb_packages/GTDB-TK_release226_2025_04_11"
THREADS=64

## Directory for GTDB-Tk
GTDBTK_DIR=$OUTPUT_DIR/GTDB-Tk
mkdir -p $GTDBTK_DIR

## Run GTDB-Tk
for assembly in "$GENOME_DIR/*.fasta"; do
    sample=$(basename "$assmebly" .fasta)
    sample_dir="$GTDBTK_DIR/$sample"
    mkdir -p "$sample_dir"

    echo "Running GTDB-Tk for $sample ....."
    
    gtdbtk classify_wf \
        --genome_file "$assembly" \
        --out_dir "$sample_dir" \
        --data_path "$GTDBTK_DATA_PATH" \
        --cpus "$THREADS"

    echo "Finished $sample"
    echo
done
