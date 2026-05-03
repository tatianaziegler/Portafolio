# Genomic Data Analysis Portfolio

## 👩‍🔬 Tatiana Ziegler, PhD

Population Genomics | Bioinformatics | Livestock Genetics

---

## 📌 Overview

This repository contains a collection of genomic data analysis projects developed using real SNP datasets from livestock populations.

The main goal of this portfolio is to demonstrate practical experience in:
- Population structure analysis
- Genetic differentiation and selection signatures
- Runs of Homozygosity (ROH) and inbreeding estimation
- Data visualization and reproducible genomic workflows

All analyses are implemented in R and based on standard bioinformatics tools used in population and animal genetics.

---

## 🧬 Projects included

### 1. FST Analysis (Selection Signatures)
Folder: `fst_analysis/`

- Calculation of genome-wide FST between POP1 and POP2
- SNP-level differentiation analysis
- Manhattan plot of selection signals
- Identification of highly differentiated loci

---

### 2. Population Structure (PCA)
Folder: `pca_analysis/`

- Principal Component Analysis using SNP genotype data
- 2D and 3D visualization of genetic structure
- Comparison between POP1 and POP2
- SNP filtering and genotype preprocessing

---

### 3. Runs of Homozygosity (ROH)
Folder: `roh_analysis/`

- Detection of ROH using sliding window approaches
- Classification of ROH by length (recent vs ancient inbreeding)
- Estimation of genomic inbreeding coefficient (FROH)
- Population comparison of inbreeding levels

---

### 4. Visualization
Folder: `visualization/`

- Manhattan plot (FST selection scan)

---

## 🧪 Data and Methods

- SNP datasets generated using PLINK workflows
- R packages used:
  - ggplot2
  - tidyverse
  - detectRUNS
  - vcfR
  - adegenet
  - plotly


---

## 📌 Notes

This portfolio is intended as a demonstration of applied bioinformatics skills for research and industry positions in:

- Animal breeding and genetics
- Population genomics
- Bioinformatics analysis pipelines

---

## 📬 Contact

Tatiana Ziegler, PhD  
Email: zieglertatiana@gmail.com  
Location: Germany / Argentina

---

## ⚠️ Reproducibility

All scripts are written for transparency and reproducibility.  
Input data formats are standard PLINK outputs (.ped, .map, .bed, .bim, .fam, .vcf).

