# Genome Assembly and Analysis Pipeline
This repository contains the data processing workflow developed as part of our 7th semester project, which aims to characterize soil-derived bacterial isolates using Oxford Nanopore sequencing and explore their potential for antibiotic production.

The pipeline is designed to process raw sequencing data, perform quality assessments, assemble genomes, polish assemblies, and evaluate final assembly quality.

---

## Repository Contents
The repository currently contains the following scripts and directories:
- `0. Full Script`: A comprehensive pipeline script that automates the entire analysis workflow from raw Nanopore FASTQ files through to assembled and polished genomes.

The following modular scripts allow you to run each data processing step independently:
- `1. Basic Statistics`: Generate basic sequencing statistics such as read counts, total bases, and average read length of each sample.
- `2. NanoPlot`: Visualise sequence run quality using NanoPlot. 
- `3. Filtering`: Filter raw reads based on quality and length thresholds to improve assembly quality.
- `4. Assembly`: Assemble filtered reads into draft genomes using Flye optimised for Nanopore data.
- `5. Polishing`: Improve assembly accuracy by iterative error correction with long reads.
- `6. Quality of Assembly`: Evaluate assembly quality metrics using QUAST.

As the project develops, more files and folders may be added.

---

## Workflow Overview
The pipeline performs the following steps:
1. **Basic statistics on raw FASTQ reads**  
   - Counting reads, bases, and calculating average read length.

2. **Quality assessment using NanoPlot**

4. **Read filtering using Filtlong**  
   - Minimum read length: 1000 bp  
   - Keep top 90% highest-quality reads

5. **Genome assembly with Flye**  
   - Using filtered Nanopore reads  
   - Genome size set to 5 Mb

6. **Polishing assemblies using Racon**  
   - Mapping reads with minimap2  
   - Sorting/indexing with Samtools  
   - Racon polishing

7. **Assembly quality assessment with QUAST**  

---

## Setup and How to Run the Pipeline
Download the script to be run. 
When a scripts is run for the first time, run this code in the Terminal
```bash
chmod +x file_name.sh
````

Place your FASTQ files in:
Data

Then run:
```bash
file_name.sh
```
### Conda Environments
Conda environments allow you to create isolated spaces with specific packages and versions, avoiding conflicts between projects.

If you have not done it yet, you need to create the needed environments: 
*Nanoplot*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*Filtlong*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*Flye*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*minimap2*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*samtools*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*racon*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*QUAST*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```

### Conda Environment Activation Issue
When running scripts inside VS Code or Code-Server using the "Run Active File" feature, you might encounter an error like:
```
CondaError: Run 'conda init' before 'conda activate'
```
This happens because VS Code runs script in a non-interactive shell. 

To fix it, add:
```bash
CONDA_BASE=$(dirname $(dirname $(which conda)))
source "$CONDA_BASE/etc/profile.d/conda.sh"
```
Before activating conda environment

## Directories


---

## Software Requirements
This pipeline relies on several bioinformatics tools and libraries. We recommend using Conda to manage dependencies easily:
- NanoPlot (v)
- Filtlong (v)
- Flye (v)
- minimap2 (v)
- samtools (v)
- racon (v)
- QUAST (v)


---

## Resources
Nanoplot:
- https://github.com/wdecoster/NanoPlot
  
Filtlong:
- https://github.com/rrwick/Filtlong

Flye:
- https://github.com/mikolmogorov/Flye

Racon:
- https://github.com/isovic/racon

QUAST:
- https://github.com/ablab/quast




