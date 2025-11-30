###################################
#!/bin/bash

# Directories
INPUT_DIR="Data"          
OUTPUT_DIR="Results"
MEDAKA_DIR="$OUTPUT_DIR/Medaka"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Conda environment
CHECKM="checkm_env"

###################################

# CheckM
echo "---------------------------------------------------"
echo "Running CheckM lineage workflow"
echo "---------------------------------------------------"

## Directory for CheckM
CHECKM_DIR="$OUTPUT_DIR/CheckM"
mkdir -p "$CHECKM_DIR"

## Activate CheckM environment
conda activate "$CHECKM"

## CheckM parameters
THREADS=8
CHECKM_DB="/databases/checkm_data_2015_01_16/"
export CHECKM_DB

## Make Genome Directory
GENOME_DIR="Results/CheckM/Genomes"
mkdir -p "$GENOME_DIR"

echo "Collecting genomes for CheckM ....."

for assembly in "$MEDAKA_DIR"/*; do
    sample=$(basename "$assembly")
    fasta="$assembly/consensus.fasta"

    if [ -f "$fasta" ]; then
        sample_dir="$GENOME_DIR/$sample"
        mkdir -p "$sample_dir"

        cp "$fasta" "$sample_dir/$sample.fasta"
        echo "Added $sample.fasta to $sample_dir"
    else
        echo "WARNING: $fasta not found"
    fi

    echo "Running CheckM for $sample ....."

    checkm lineage_wf \
        --reduced_tree \
        -t $THREADS \
        -x fasta \
        "$sample_dir" \
        "$CHECKM_DIR/$sample"

    echo "CheckM completed for $sample ....."

done
conda deactivate
