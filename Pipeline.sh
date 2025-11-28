###################################
#!/bin/bash

# Input and output directory
INPUT_DIR="Data"
OUTPUT_DIR="Results"

# Load conda commands
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Conda environments
## check in terminal: 'conda env list'
NANOPLOT="nanoplot2"
FILTLONG="filtlong_env"
FLYE="flye_env"
MINIMAP2="minimap2"
SAMTOOLS="samtools_env"
MEDAKA="medaka"
QUAST="quast_env"
CHECKM=
CHECKM2="checkm2_env"
BAKTA="bakta1.10.3_env"

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

###################################

# Assess read quality of filtered reads with NanoPlot
echo "---------------------------------------------------"
echo "Assessing read quality of filtered reads with NanoPlot after filtering"
echo "---------------------------------------------------"

## Load NanoPlot environment
conda activate $NANOPLOT

## Directory for NanoPlot Results
NANOPLOT_DIR="$OUTPUT_DIR/NanoPlot_filtered"
mkdir -p "$NANOPLOT_DIR"

## Run NanoPlot
for file in "$FILTERED_DIR"/*/*_filtered.fastq; do
    sample=$(basename "$file" _filtered.fastq)
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

###################################

# Assembly of sequencing data
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

###################################


# Calculate coverage for all assemblies
echo "---------------------------------------------------"
echo "Calculating coverage"
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
    conda activate $MINIMAP2
    minimap2 -a -x map-ont "$current_assembly" $reads > "$COVERAGE_DIR/$sample.bam"
    conda deactivate

    ## Sort and index
    conda activate $SAMTOOLS
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

###################################

# Polishing with Medaka
echo "---------------------------------------------------"
echo "Polishing with Medaka"
echo "---------------------------------------------------"

## Medaka parameters
THREADS=16        
MODEL="r1041_e82_400bps_hac_v5.0.0"        

## Directory for Medaka
MEDAKA_DIR="$OUTPUT_DIR/Medaka"
mkdir -p "$MEDAKA_DIR"

## Load Medaka environment
conda activate $MEDAKA

## Run Medaka on assemblies
for polished_assembly in "$ASSEMBLY_DIR"/*; do
    sample=$(basename "$polished_assembly")
    sample_dir="$MEDAKA_DIR/$sample"
    mkdir -p "$sample_dir"

    echo "Starting polishing with Medaka for $sample ....."

    polished_assembly="$ASSEMBLY_DIR/$sample/assembly.fasta"
    reads="$FILTERED_DIR/$sample/${sample}_filtered.fastq"
    
    medaka_consensus \
        -i $reads \
        -d $polished_assembly \
        -o $sample_dir \
        -t $THREADS \
        -m $MODEL

    echo
    echo "Completed Medaka for $sample. Results have been saved in "$MEDAKA_DIR""
    echo "---------------------------------------------------"
done
conda deactivate

###################################

# Run QUAST on output
echo "---------------------------------------------------"
echo "Assessing assembly quality with QUAST"
echo "---------------------------------------------------"

## Load QUAST environment
conda activate $QUAST

## Directory for QUAST
QUAST_DIR="$OUTPUT_DIR/QUAST"
mkdir -p "$QUAST_DIR"

## QUAST parameters
THREADS=8        # change

## Run QUAST
for assembly in "$ASSEMBLY_DIR"/*; do
    [ -d "$assembly" ] || continue
    sample=$(basename "$assembly")

    raw_assembly="$ASSEMBLY_DIR/$sample/assembly.fasta"
    medaka_polished="$MEDAKA_DIR/$sample/consensus.fasta" 

    echo "Running QUAST for $sample ....."
    
    if [[ -f "$medaka_polished" ]]; then
        quast.py \
            "$raw_assembly" "$medaka_polished" \
            -o "$QUAST_DIR/${sample}" \
            -l "${sample}_raw,${sample}_medaka" \
            --threads "$THREADS"
    else
        echo "Medaka polished file NOT found for $sample – running QUAST without medaka"
        quast.py \
            "$raw_assembly" \
            -o "$QUAST_DIR/${sample}" \
            -l "${sample}_raw" \
            --threads "$THREADS"
    fi
    
    echo "QUAST report saved in: $QUAST_DIR"
    echo "---------------------------------------------------"
done
conda deactivate

###################################

# Run CheckM2
echo "---------------------------------------------------"
echo "Assessing completeness and contamination on nonpolished assemblies with CheckM2"
echo "---------------------------------------------------"

## Load CheckM2 environment
conda activate $CHECKM2

## Directory for CheckM2 
CHECKM2_DIR="$OUTPUT_DIR/CheckM2"
mkdir -p "$CHECKM2_DIR"

## CheckM2 database path
CHECKM2_DB="Database/CheckM2_database/uniref100.KO.1.dmnd"
export CHECKM2DB="$CHECKM2_DB"

## CheckM2 parameters
THREADS=16      #change
ASSEMBLIES="$ASSEMBLY_DIR" #change if needed
echo "Running CheckM2 on nonpolished assemblies in $ASSEMBLIES"

## Run CheckM2
for fasta_file in "$ASSEMBLIES"/*/*.fasta;  do
    sample=$(basename "$(dirname "$fasta_file")")
    sample_dir="$CHECKM2_DIR/$sample/Nonpolished"
    mkdir -p "$sample_dir"

    echo "Running CheckM2 on sample: $sample"
    checkm2 predict \
        --threads $THREADS \
        --input "$fasta_file" \
        --output_directory "$sample_dir" \
        --database_path "$CHECKM2_DB"

    echo "Results saved in $sample_dir"
    echo "----------------------------"
done
conda deactivate

###################################

# Run CheckM2
echo "---------------------------------------------------"
echo "Assessing completeness and contamination of polished assemblies with CheckM2"
echo "---------------------------------------------------"

## Load CheckM2 environment
conda activate $CHECKM2

## Directory for CheckM2 
CHECKM2_DIR="$OUTPUT_DIR/CheckM2_medaka"
mkdir -p "$CHECKM2_DIR"

## CheckM2 database path
CHECKM2_DB="Database/CheckM2_database/uniref100.KO.1.dmnd"
export CHECKM2DB="$CHECKM2_DB"

## CheckM2 parameters
THREADS=16      #change
ASSEMBLIES="$MEDAKA_DIR" #change if needed
echo "Running CheckM2 on polished assemblies in $ASSEMBLIES"

## Run CheckM2
for fasta_file in "$ASSEMBLIES"/*/*.fasta;  do
    sample=$(basename "$(dirname "$fasta_file")")
    sample_dir="$CHECKM2_DIR/$sample/Medaka"
    mkdir -p "$sample_dir"

    echo "Running CheckM2 on sample: $sample"
    checkm2 predict \
        --threads $THREADS \
        --input "$fasta_file" \
        --output_directory "$sample_dir" \
        --database_path "$CHECKM2_DB"

    echo "Results saved in $sample_dir"
    echo "----------------------------"
done
conda deactivate

###################################



