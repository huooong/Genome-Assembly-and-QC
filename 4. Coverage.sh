###################################
#!/bin/bash

# Directories
INPUT_DIR="Data"
OUTPUT_DIR="Results"
ASSEMBLY_DIR="$OUTPUT_DIR/Assembly"
FILTERED_DIR="$OUTPUT_DIR/Filtered"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

###################################

# Calculate coverage for all assemblies
echo "---------------------------------------------------"
echo "Calculating coverage for all samples"
echo "---------------------------------------------------"

## Directory for coverage analysis
COVERAGE_DIR="$OUTPUT_DIR/Coverage"
mkdir -p "$COVERAGE_DIR"

## Output file
COVERAGE_FILE="$COVERAGE_DIR/coverage_summary.tsv"
echo -e "Sample\tMean_coverage" > "$COVERAGE_FILE"

## Run coverage analysis for all assemblies
for assembly in "$ASSEMBLY_DIR"/*; do
    sample=$(basename "$assembly")
    current_assembly="$ASSEMBLY_DIR/$sample/assembly.fasta"
    reads="$FILTERED_DIR/$sample/"*_filtered.fastq

    echo "Processing $sample ....."
   
    ## Map reads to assembly
    conda activate minimap2
    minimap2 -a -x map-ont "$current_assembly" $reads > "$COVERAGE_DIR/$sample.bam"
    conda deactivate

    ## Sort and index
    conda activate samtools_env
    samtools sort -o "$COVERAGE_DIR/$sample.sorted.bam" "$COVERAGE_DIR/$sample.bam"
    samtools index "$COVERAGE_DIR/$sample.sorted.bam"

    ## Coverage depth file
    samtools depth -a "$COVERAGE_DIR/$sample.sorted.bam" > "$COVERAGE_DIR/$sample.depth.txt"
    
    ## Compute mean coverage
    mean_coverage=$(awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}' "$COVERAGE_DIR/$sample.depth.txt")

    ## Write to summary
    echo -e "$sample\t$mean_coverage" >> "$COVERAGE_FILE"

    echo "Mean coverage for $sample: $mean_coverage"
    conda deactivate
done

# Remove intermediate files
rm "$COVERAGE_DIR/"*.bam
rm "$COVERAGE_DIR/"*.bam.bai
rm "$COVERAGE_DIR/"*.txt


