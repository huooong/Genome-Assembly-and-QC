# Soil Bacteria Genome Analysis – 7th Semester Project

This repository contains the data processing workflow used in a 7th semester biotechnology project focused on characterizing soil-derived bacterial isolates and exploring their potential for antibiotic production.  
The analyses include read quality assessment, filtering, genome assembly, polishing, and assembly evaluation.

---

## Repository Contents

Currently the repository includes:

- `P7_script.sh` – A full analysis pipeline from raw Nanopore reads to assembled and polished genomes.

As the project develops, the following folders may be added:
/data/ # Raw and processed sequencing data
/results/ # NanoPlot, assemblies, polished genomes, QUAST reports
/scripts/ # Modular shell scripts (if workflow is split)


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

Place your FASTQ files in:
P7/Data/


Then run:
```bash
bash P7_script.sh
```

Results will be saved in:
P7/Results/

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




