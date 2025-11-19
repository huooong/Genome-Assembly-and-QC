# Genome Assembly and Analysis of Soil Bacteria

This repository contains the data processing workflow used in our 7th semester project titled 'Identification of Next Generation Antibiotics using Genome Mining Tools', focused on characterizing soil-derived bacterial isolates and exploring their potential for antibiotic production.  

---

## Repository Contents

Currently the repository includes:
- `0. Full Script` – A combined full analysis pipeline from raw Nanopore reads to assembled and polished genomes.

The scripts listed below are structured to allow each processing step to be run separately:
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

### Conda Environments
Conda environments allow you to create isolated spaces with specific packages and versions, avoiding conflicts between projects.

To create the needed environments: 


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




