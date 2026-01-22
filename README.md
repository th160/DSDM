# DSDM

[![R](https://img.shields.io/badge/Made_with-R-blue.svg)](https://cran.r-project.org/)

Enabling DNA Intrinsic Topology  for Genome Annotations by DSDM

## Overview

DSDM (DNA Shape Dependent Mixture) is a genome annotation model based on intrinsic DNA topology. Firstly, the input sequencing data, including the DNase-seq data and the ATAC-seq data, are collected to extract highly ranked peak regions present open chromatin regions which are accessible to transcription factors and other regulatory proteins. Then, the DNA shape features of extracted peak summits are profiled at the nucleotide level. Next, we constructed the dependent mixture model based on extracted shape features. After training, the genome annotation work can be discovered using our model, especially the motif discovery task. Simultaneously, through the DSDM model, we can obtain the local DNA topological landscape for given DNA sequences. This repository contains the source code for data preparation, model training, and application of the proposed model.

### Workflow Visualization

![Model Workflow Diagram](workflow.png)
---


## 🛠️ Prerequisites
This project is built using R. Please ensure your environment meets the following requirements:
* **R Version**: >= 4.0.0
### Required Packages
```r
# Install required CRAN packages
install.packages("depmixS4")
library(depmixS4)
```

### Project Structure
```r
.
├── README.md             
├── dataset/              # Folder containing raw and processed data files
├── model/                # Folder for pre-trained model parameters
└── scripts/              # Source codes for training and application
    ├── DSDM_modeling.R              # The main script used for training the model
└── DSDM_prediction/ 
    ├── DSDM_prediction.R             # Application script to apply the trained model to new data
```

Due to GitHub's file size limitations, the whole pre-trained model are too large to be hosted directly in this repository. Pre-trained DSDM Model and processed data via **[Zenodo](https://doi.org/10.5281/zenodo.17896569)**
