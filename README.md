# Genome Analysis Pipeline
This repository contains the data processing workflow developed as part of our 7th semester project, which aims to characterize soil-derived bacterial isolates using Oxford Nanopore sequencing and explore their potential for antibiotic production.

The pipeline is designed to process raw sequencing data, perform quality assessments, assemble genomes, polish assemblies, and evaluate final assembly quality.

---

## Repository Contents
The repository currently contains the following scripts and directories:
- `0. Full Script`: (In progress) A comprehensive pipeline script that automates the entire analysis workflow from raw Nanopore FASTQ files through to assembled and polished genomes.

The following modular scripts allow you to run each data processing step independently:
- `1. Basic Statistics`
- `2.1. NanoPlot on raw reads`
- `2.2. NanoPlot on filtered reads`
- `3. Filtering`
- `4. Assembly`
- `5.1. Polishing with Racon`
- `5.2. Polishing with Medaka`
- `6.1. QUAST`
- `6.2. BUSCO` (In progress)
- `7. Completeness and contamination`
- `8. Annotation` (In progress)
- `MultiQT`

As the project develops, more files and folders may be added.

---

## Workflow Overview
The pipeline performs the following steps:

**1.0. Basic statistics on raw FASTQ reads**  
   - Generate basic sequencing statistics for each samples, including:
      - Read counts,
      - Total bases,
      - Average read length.

**2.1. Quality assessment of raw reads using NanoPlot**
   - Provides summary statistics and visualisations describing the quality of unfiltered reads.
   - Generated plots include:
      - Weighted and non-weighted histogram of read lengths,
      - Length vs. Quality scatter plot,
      - Yield by length plot.

**2.2. Quality assessment of filtered reads using NanoPlot**
   - Runs NanoPlot on Filtlong-filtered reads to evaluate quality after filtering.

**3.0. Read filtering using Filtlong**  
   - Filters raw Nanopore reads to improve downstream assembly quality using the following parameters:
      - Minimum read length: **1000 bp**,  
      - Keep top **90%** highest-quality reads,
      - Trim bases until only the best **500 Mbp** remain.

**4.0. Genome assembly with Flye**  
   - Assembles genomes from filtered Nanopore reads using the following parameters:
      - Estimated genome size: **5 Mb**.
      - Threads: **4**.

**5.1 Polishing assemblies using Racon**  
   - Performs iterative polishing (**3 iterations**) consisting of:
      - Mapping reads with minimap2.  
      - Polishing with Racon.

**5.2 Final polish using Medaka**  
   - Refines assemblies using Medaka with the model:
      - `r1041_e82_400bps_hac_v5.0.0`

**6.0. Assembly quality assessment with QUAST**
   - Evaluates assembly statistics including N50, genome size, misassemblies, and GC content.
  
**7.0. Assessment of completeness and contamination with CheckM2**
   - Estimates genome completeness and contamination to assess assembly quality and reliability.

**8.0. MultiQT**
   - Combines key results and quality metrics from multiple tools into a single, interactive HTML report.

---
## Directories 
### Project Structure
This project is organised to streamline the processing and analysis of sequencing data. 

The main project directory is called `P7`. Within `P7`, there is a `Data` directory containing all raw FASTQ files (`*.fastq.gz`) generated from Oxford Nanopore sequencing. 

The fastq.gz files are named `PBI54877.barcodeX.fastq.gz`, where X = 17-36, corresponding to 20 individual samples.

### Output Organisation
All analysis results are stored inside a folder called `Results`. This directory contains subfolders for each step of the pipeline: 
   - Filtered/
   - NanoPlot/
   - NanoPlot_filtered/
   - Assembly/
   - Racon/
   - Medaka/
   - QUAST/
   - CheckM2/
   - MultiQT/

Each of these subfolders contains one directory per sample, following the same naming scheme as the input files. 

---

## Setup and How to Run the Pipeline
1. Download/copy the script.
   
2. Make the script executable.
Before running a script for the first time, execute the following code:
```bash
chmod +x file_name.sh
````
Replace `file_name.sh` with the specific script. 

3. Prepare input files.
Place all raw FASTQ files in the directory `Data`.
```bash
cp -r /raw_data/MA/tinyearth/20251119_np_PBI54877/*.fastq.gz ~/P7/Data/
``` 

5. Run the script.
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
*Medaka*
```bash
conda create -n medaka -c conda-forge -c nanoporetech -c bioconda medaka
conda install -c conda-forge pyabpoa
```
*QUAST*
```bash
conda create -n quast_env -c bioconda quast
```
*CheckM2*
```bash
conda create -n checkm2_env -c bioconda -c conda-forge checkm2
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

### Software Requirements
This pipeline relies on several bioinformatics tools and libraries. We recommend using Conda to manage dependencies easily:
- NanoPlot (v1.46.1)
- Filtlong (v0.2.1)
- Flye (v2.9.6-b1802)
- minimap2 (v2.30-r1287)
- samtools (v1.22.1)
- racon (v1.5.0)
- QUAST (v5.3.0)
- Medaka (v2.1.1)
- CheckM2 (v1.1.0)

---
## CheckM2
In order to use CheckM2 the DIAMOND database needs to be downloaded. This is done by:
```bash
checkm2 database --download
```
I moved the database folder to be inside my working directory, but you can do what you want to, just make sure to put down the right path in the script.

Testrun can be done using:
```bash
checkm2 testrun
```

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

Medaka: 
- https://github.com/nanoporetech/medaka

CheckM2:
- https://github.com/chklovski/CheckM2/




