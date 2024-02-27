
# README for Key Pathway Miner (KPM) Analysis Pipeline

## Overview
This repository contains a Nextflow pipeline for Key Pathway Miner (KPM) analysis, aimed at identifying key pathways in proteomics data. The pipeline integrates differential expression data with network information to pinpoint significant pathways under different conditions and timepoints.

## Usage
To run the pipeline, specify the following flags:

1. `--meta_file`: Path to the metadata file, detailing sample information like experimental conditions.
2. `--count_file`: Path to the count file, containing expression data for analysis.
3. `--network_file`: Path to the network file, which provides the interaction network information for pathway analysis.
4. `--logFC`: Boolean flag to apply a log fold change (logFC) threshold (default: `true`).
5. `--logFC_up`: Upper log2 fold change threshold for determining upregulated elements (default: `1`).
6. `--logFC_down`: Lower log2 fold change threshold for downregulated elements (default: `-1`).
7. `--p_adj`: Boolean flag to use adjusted p-values (default: `true`).
8. `--alpha`: Significance threshold for p-values (default: `0.05`).
9. `--output`: Directory for storing output files (default: `"./output/"`).

Example command:
```
./nextflow run main.nf --meta_file path/to/meta.txt --count_file path/to/count.txt --network_file path/to/network.txt --output output_directory
```

## Output Description
The KPM analysis pipeline generates various outputs, including:

- **Pathway Analysis Results**: Detailed results of the pathway analysis, including key pathways identified under different conditions.
- **Statistical Summaries**: Summaries of the statistical analysis, including p-values and fold change information for elements in the identified pathways.
- **Visualization Files**: Graphical representations of the key pathways, aiding in the interpretation and presentation of the results.

## Getting Started
Begin by cloning this repository and ensuring Nextflow is installed on your system. Prepare your metadata, count, and network files according to the required format. Execute the command with the correct paths to your input files. The results will be available in the specified output directory.

---
**Note**: Adjust paths, filenames, and parameters as needed to fit your project's specific requirements and data structure.
