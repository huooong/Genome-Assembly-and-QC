###################################
#!/bin/bash

# Directories
INPUT_DIR="Results"
OUTPUT_DIR="Results"
MEDAKA_DIR="$OUTPUT_DIR/Medaka"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Conda environment
BAKTA="bakta_env"

###################################

# Genome annotation with BAKTA
echo "---------------------------------------------------"
echo "Running BAKTA for gene annotation"
echo "---------------------------------------------------"

## Directory for BAKTA 
BAKTA_DIR="$OUTPUT_DIR/Bakta"
mkdir -p "$BAKTA_DIR"

## Activate bakta environment
conda activate $BAKTA

## Bakta parameters
THREADS=8
BAKTA_DB="databases/bakta/20240125"

## Run Bakta on all assemblies
for polished_assembly in "$MEDAKA_DIR"/*; do
    sample=$(basename "$polished_assembly")
    sample_dir="$BAKTA_DIR/$sample"
    mkdir -p "$sample_dir"
    medaka_polished="$MEDAKA_DIR/$sample/consensus.fasta"

    echo "Annotating $sample ....."

    bakta \
      --db "$BAKTA_DB \
      --threads $THREADS \
      --output "$sample_dir" \
      "$medaka_polished"

    echo "Done with $sample"
    echo "---------------------------------------------------"

done
conda deactivate



