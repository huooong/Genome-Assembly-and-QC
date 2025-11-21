# Genome Analysis Pipeline
This repository contains the data processing workflow developed as part of our 7th semester project, which aims to characterize soil-derived bacterial isolates using Oxford Nanopore sequencing and explore their potential for antibiotic production.

The pipeline is designed to process raw sequencing data, perform quality assessments, assemble genomes, polish assemblies, and evaluate final assembly quality.


---


## Repository Contents
The repository currently contains the following scripts and directories:
- `0. Full Script`: A comprehensive pipeline script that automates the entire analysis workflow from raw Nanopore FASTQ files through to assembled and polished genomes.

The following modular scripts allow you to run each data processing step independently:
- `1. Basic Statistics`
- `2. NanoPlot`
- `3. Filtering`
- `4. Assembly`
- `5. Polishing`
- `6. Quality of Assembly`

As the project develops, more files and folders may be added.


---


## Workflow Overview
The pipeline performs the following steps:
1. **Basic statistics on raw FASTQ reads**  
   - Generate basic sequencing statistics such as read counts, total bases, and average read length of each sample.

2. **Quality assessment using NanoPlot**
   - Gives summary statistics and plots describing the quality of input data prior to analysis.
   - Plots that are generated, include:
      - Weighted and non-weighted histogram of read lengths.
      - Length vs. Quality scatter plot.
      - Yield by length plot.

4. **Read filtering using Filtlong**  
   - Filter raw reads based on quality and length thresholds to improve assembly quality.
      - Minimum read length: 1000 bp.  
      - Keep top 90% highest-quality reads.

5. **Genome assembly with Flye**  
   - Using filtered Nanopore reads.  
   - Genome size set to 5 Mb.
   - Threads set to 4.

6. **Polishing assemblies using Racon**  
   - Mapping reads with minimap2.  
   - Sorting/indexing with Samtools.  
   - Racon polishing.

7. **Assembly quality assessment with QUAST**
   
8. **MultiQT**


---


## Setup and How to Run the Pipeline

1. Download the script.
   
3. Make the script executable.
Before running a script for the first time, execute the following code:
```bash
chmod +x file_name.sh
````
Replace `file_name.sh` with the specific script. 

3. Prepare input files.
Place all raw FASTQ files in the directory `Data`.
Ensure the files end with `.fastq`.

4. Run the script.
```bash
./file_name.sh
```
The script will process all FASTQ files found in the `Data`-directory and output results to automatically created result folders.

---

### Conda Environments
Conda environments allow you to create isolated spaces with specific packages and versions, avoiding conflicts between projects.

If you have not done it yet, you need to create the needed environments: 
*Nanoplot*
```bash
conda create -n nanoplot2 -c bioconda nanoplot
```
*Filtlong*
```bash
conda create -n filtlong_env -c bioconda filtlong
```
*Flye*
```bash
conda create -n flye_env -c bioconda flye
```
*minimap2*
```bash
conda create -n minimap2 -c bioconda minimap2
```
*samtools*
```bash
conda create -n samtools_env -c bioconda samtools
```
*racon*
```bash
conda create -n racon_env -c bioconda racon=1.4.2
```
*QUAST*
```bash
conda create -n quast_env -c bioconda quast
```
*MultiQT*
```bash
conda create -n multiqc -c bioconda multiqc
```

---

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

---

## Directories


---

## Software Requirements
This pipeline relies on several bioinformatics tools and libraries. We recommend using Conda to manage dependencies easily:
- NanoPlot (v1.46.1)
- Filtlong (v0.2.1)
- Flye (v2.9.6-b1802)
- minimap2 (v2.30-r1287)
- samtools (v1.22.1)
- racon (v1.5.0)
- QUAST (v5.3.0)

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




