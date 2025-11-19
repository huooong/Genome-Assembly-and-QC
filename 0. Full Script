#!/bin/bash

#Input and output directory (CORRECT THIS IF NEEDED)
INPUT_DIR="P7/Data"
OUTPUT_DIR="P7/Results"

#Basic Statistics on the Sequencing Data
echo "STEP 1: Counting reads, bases and average read length"

for file in "$INPUT_DIR"/*.fastq; do
    sample=$(basename "$file" .fastq)
    echo "Processing $sample..."

    ## Count reads
    reads=$(awk '{s++} END{print s/4}' "$file")

    ## Count bases
    bases=$(awk 'NR%4==2 {bases+=length($0)} END{print bases}' "$file")

    ## Average read length
    avg_length=$(awk 'NR%4==2 {bases+=length($0); reads++} END{print bases/reads}' "$file")

    ## Print results
    echo " Reads: $reads"
    echo " Bases: $bases"
    echo " Average read length: $avg_length"
    echo "-------------------------------------"
done

# Run NanoPlot
echo "STEP 2: Assessing quality with NanoPlot"

##Load NanoPlot environment
conda activate nanoplot2

##Create directory for NanoPlot Results
NANOPLOT_OUT="$OUTPUT_DIR/NanoPlot"
mkdir -p "$NANOPLOT_OUT"

for file in "$INPUT_DIR"/*.fastq; do
    sample=$(basename "$file" .fastq)
    echo "Running NanoPlot on $sample..."
    NanoPlot --fastq "$file" -o NanoPlot --prefix "$sample" --tsv_stats --raw
done
conda deactivate
echo "NanoPlot results have been saved in $NANOPLOT_DIR"

# Filter reads using Filtlong
echo "STEP 3: Filtering reads using Filtlong"

##Load Filtlong environment
conda activate filtlong_env

##Directory for filtered results
FILTERED_DIR="$OUTPUT_DIR/Filtered"
mkdir -p "$FILTERED_DIR"

for file in "$INPUT_DIR"/*.fastq; do
    sample=$(basename "$file" .fastq)
    echo "Filtering $file ..."

    filtlong --min_length 1000 --keep_percent 90 "$file" > "${FILTERED_DIR}/${sample}_filtered.fastq"
done
conda deactivate
echo "All files have been filtered. Filtered reads have been saved in $FILTERED_DIR"

# Assembling of sequencing data
echo "STEP 4: Assembling sequenced data using Flye Assembler"

##Load Flye environment
conda activate flye_env

##Create directory for assembly
ASSEMBLY_DIR="$OUTPUT_DIR/Assembly"
mkdir -p "ASSEMBLY_DIR"

##Flye parameters
GENOME_SIZE="5m"
THREADS=4

##Run flye
for file in "$FILTERED_DIR"/*_filtered.fastq; do
	sample=$(basename "$file"_filtered.fastq)
	echo "Assembling $sample..."
flye --nano-raw "$file" --genome-size "$GENOME_SIZE" --out-dir "ASSEMBLY_DIR/$sample" --threads "$THREADS"
done
conda deactivate
echo “Assemblies saved in $ASSEMBLY_DIR”

# Polishing using Racon
echo "STEP 5: Polishing assemblies with Racon"

##Create directory for polished assemblies
POLISH_DIR="$OUTPUT_DIR/Polished"
mkdir -p "$POLISH_DIR"

for assembly_folder in "$ASSEMBLY_DIR"/*; do
    sample=$(basename "$assembly_folder")
    assembly="$assembly_folder/assembly.fasta"
    reads="$FILTERED_DIR/${sample}_filtered.fastq"

    echo "Polishing $sample ..."

    ##Map reads back to the assembly
    conda activate minimap2
    minimap2 -x map-ont -t "$THREADS" "$assembly" "$reads" > "$POLISH_DIR/${sample}_aln.sam"
    conda deactivate
    

    ##Sort and index alignment
    conda activate samtools_env
    samtools sort "$POLISH_DIR/${sample}_aln.sam" -o "$POLISH_DIR/${sample}_aln_sorted.bam"
    samtools index "$POLISH_DIR/${sample}_aln_sorted.bam"
    conda deactivate
    

    ##Polish with Racon
    conda activate racon_env
    racon -t "$THREADS" "$reads" "$POLISH_DIR/${sample}_aln.sam" "$assembly" > "$POLISH_DIR/${sample}_racon.fasta"
    conda deactivate
    
    echo "$sample polished assembly saved at $POLISH_DIR/${sample}_racon.fasta"
done

# Run QUAST on output
echo "STEP 6: Assessing assembly quality with QUAST"
conda activate quast_env

##Create directory for QUAST
QUAST_OUT="$OUTPUT_DIR/QUAST"
mkdir -p "$QUAST_OUT"

#Run QUAST 
for folder in "$ASSEMBLY_DIR"/*; do
    [ -d "$folder" ] || continue
    sample=$(basename "$folder")

    raw_assembly="$ASSEMBLY_DIR/$sample/assembly.fasta"
    polished_assembly="$OUTPUT_DIR/Polished/${sample}_racon.fasta"

    echo "Running QUAST for $sample ..."
    quast.py \
        "$raw_assembly" "$polished_assembly" \
        -o "$QUAST_OUT/${sample}" \
        -l "${sample}_raw,${sample}_racon" \
        --threads "$THREADS"
done

conda deactivate
echo "QUAST reports saved in: $QUAST_OUT"
