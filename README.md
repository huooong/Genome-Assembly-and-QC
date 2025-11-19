# Genome Assembly and Analysis of Soil Bacteria

This repository contains the data processing workflow used in our 7th semester project titled 'Identification of Next Generation Antibiotics with Genome Mining', focused on characterizing soil-derived bacterial isolates and exploring their potential for antibiotic production.  

---

## Repository Contents

Currently the repository includes:
- `0. Full Script` – A combined full analysis pipeline from raw Nanopore reads to assembled and polished genomes.

The scripts listed below are structured to allow each processing step to be run separately:
- `1. Basic Statistics`
- `2. NanoPlot`

As the project develops, more files and folders may be added.
---

## Workflow Overview

The pipeline performs the following steps:
1. **Basic statistics on raw FASTQ reads**  
   - Counting reads, bases, and calculating average read length.

2. **Quality assessment using NanoPlot**

3. **Read filtering using Filtlong**  
   - Minimum read length: 1000 bp  
   - Keep top 90% highest-quality reads

4. **Genome assembly with Flye**  
   - Using filtered Nanopore reads  
   - Genome size set to 5 Mb

5. **Polishing assemblies using Racon**  
   - Mapping reads with minimap2  
   - Sorting/indexing with Samtools  
   - Racon polishing

6. **Assembly quality assessment with QUAST**  
   - Comparing raw Flye assembly and Racon-polished assembly

---

## How to Run the Pipeline
When a file is run for the first time, run this code in the Terminal.
```bash
chmod +x file.sh
````

Make sure to start in the P7 directory
```bash
cd P7
```

Place your FASTQ files in:
P7/Data/

Then run:
```bash
file.sh
```

Results will be saved in:
P7/Results/*

---

## Software Requirements
This pipeline expects the following tools to be installed, ideally through conda environments:
- NanoPlot
- Filtlong
- Flye
- minimap2
- samtools
- racon
- QUAST

---

## Conda



