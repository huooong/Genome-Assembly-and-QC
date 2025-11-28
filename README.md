# Genome Analysis Pipeline
This repository contains the data processing workflow developed as part of our 7th semester project, which aims to characterize soil-derived bacterial isolates using Oxford Nanopore sequencing and explore their potential for antibiotic production.

The pipeline is designed to process raw sequencing data, perform quality assessments, assemble genomes, polish assemblies, and evaluate final assembly quality.

---

## Repository Contents
The repository currently contains the following scripts:
- `Pipeline.sh`: A pipeline script that automates the analysis workflow from raw Nanopore FASTQ files through to assembled and polished genomes.
- `Pipeline.sbatch`: Allows you to submit the pipeline as a SLURM job.

The following modular scripts allow you to run each data processing step independently:
- `1. Filtering`
- `2.1. NanoPlot on raw reads`
- `2.2. NanoPlot on filtered reads`
- `3. Assembly with Flye`
- `4. Coverage`
- `5. Polishing with Medaka`
- `6. QUAST`
- `7. CheckM`
- `8. CheckM2`
- `Annotation` (In progress)
- `MultiQC` 
- `GTDB-Tk` (In pregress)

---

## Workflow Overview
The pipeline performs the following steps:

**Read filtering using Filtlong**  
   - Filters raw Nanopore reads to improve downstream assembly quality using the following parameters:
      - Minimum read length: **1000 bp**,  
      - Keep top **90%** highest-quality reads,
      - Trim bases until only the best **500 Mbp** remain.

**Quality assessment of raw reads using NanoPlot**
   - Provides summary statistics and visualisations describing the quality of unfiltered and filtered reads.
   - Generated plots include:
      - Weighted and non-weighted histogram of read lengths,
      - Length vs. Quality scatter plot,
      - Yield by length plot.
      - 
**Genome assembly with Flye**
   - Assembles genomes from filtered Nanopore reads using the following parameters:
      - Estimated genome size: **5 Mb**.
      - Threads: **4**.

**Calculating Coverage**  
   - Sequencing coverage for assembled genomes are calculated using minimap2 and samtools.  

**Polish assembly using Medaka**  
   - Refines assemblies using Medaka with the model:
      - `r1041_e82_400bps_hac_v5.0.0`

**Assembly quality assessment with QUAST**
   - Evaluates assembly statistics including N50, genome size, misassemblies, and GC content.
  
**Assessment of completeness and contamination with CheckM and CheckM2**
   - Estimates genome completeness and contamination to assess assembly quality and reliability.

**MultiQC**
   - Combines key results and quality metrics from multiple tools into a single, interactive HTML report.

---
## Directories 
### Project Structure
This project is organised to streamline the processing and analysis of sequencing data. 

The main project directory is called `Working_Directory`. Within `Working_Directory`, there is a `Data` directory containing all raw FASTQ files (`*.fastq.gz`) generated from Oxford Nanopore sequencing. 

The fastq.gz files are named `PBI54877.barcodeX.fastq.gz`, where X = 17-36, corresponding to 20 individual samples.

### Output Organisation
All analysis results are stored inside a folder called `Results`. This directory contains subfolders for each step of the pipeline: 
   - Filtered/
   - NanoPlot/
   - NanoPlot_filtered/
   - Assembly/
   - Coverage/
   - Medaka/
   - QUAST/
   - CheckM/
   - CheckM2/
   - MultiQC/

Each of these subfolders contains one directory per sample, following the same naming scheme as the input files. 

---

## Setup and How to Run the Pipeline
1. Download/copy the `Pipeline.sh` and `Pipeline.sbatch` to `Working_Directory`


2. Prepare input files.
Place all raw FASTQ files in the directory `Working_Directory/Data`.
```bash
cp -r /raw_data/MA/tinyearth/20251119_np_PBI54877/*.fastq.gz ~/Working_Directory/Data/
``` 
Make sure that only the samples you want to analyse are in the `Data` directory. 

3. Go to the terminal, log on, and change to your working directory:
```bash
cd Working_Directory
````

4. Run the pipeline:
```bash
sbatch Pipeline.sbatch
````

5. Status:
```bash
squeue -u $USER
```

```bash
tail -f logs/genome_assembly_and_QC_<jobID>.out
```

---

## Conda Environments
Conda environments allow you to create isolated spaces with specific packages and versions, avoiding conflicts between projects.

You can check your existing environments using: 
```bash
conda env list
```

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

*BUSCO*
```bash
conda create -n busco_env -c conda-forge -c bioconda busco=6.0.0
```

*CheckM2*
```bash
conda create -n checkm2_env -c bioconda -c conda-forge checkm2
```

*Bakta*
```bash
conda create -n bakta1.10.3_env -c conda-forge -c bioconda bakta=1.10.3
```

*MultiQC*
```bash
conda create -n multiqc -c bioconda multiqc
```

*GTDB-Tk*
```bash
conda create -n gtdbtk_env -c conda-forge -c bioconda gtdbtk=2.5.2
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

## Software Requirements
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
- Bakta (v1.10.3)
- GTDB-Tk (v2.5.2)

  
---

## Resources
Nanoplot:
- https://github.com/wdecoster/NanoPlot
  
Filtlong:
- https://github.com/rrwick/Filtlong 

Flye:
- https://github.com/mikolmogorov/Flye
- Articles:
   - **Single genome assembly:** Kolmogorov et al. (2019). Assembly of long, error-prone reads using repeat graphs. _Nature Biotechnology_ 37 (2019) 540-546.
   - **Core algorithm:** Lin et al. (2016). Assembly of long error-prone reads using de Bruijn graphs. _PNAS_ 113 (52) E8396-E8405.

Racon:
- https://github.com/isovic/racon

QUAST:
- https://github.com/ablab/quast

BUSCO:
- https://busco.ezlab.org/

Medaka: 
- https://github.com/nanoporetech/medaka

CheckM2:
- https://github.com/chklovski/CheckM2/


---
## If you want to play around with the modular scripts in VS Code: 
You can make the scripts executable. This is done by executing the following code before running a script for the first time:
```bash
chmod +x file_name.sh
````
Replace `file_name.sh` with the specific script. 

To run the script.
```bash
./file_name.sh
```



