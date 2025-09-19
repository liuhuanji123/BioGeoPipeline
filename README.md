# 🪲 Phylogenetic and Biogeographic Analysis Pipeline for Coleopteran Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a semi-automated pipeline designed for conducting phylogenetic and biogeographic analyses on families within the order **Coleoptera** (beetles).

The core methodology combines the robust phylogenetic signal from **mitogenomes** with the extensive species coverage of **DNA barcodes**. This hybrid approach produces a comprehensive and reliable **time-calibrated phylogeny**, which serves as the foundation for all subsequent biogeographic reconstructions.

The workflow integrates several key bioinformatics tools, including **RAxML-NG**, **treePL**, and **PASTML** (executed via Windows Subsystem for Linux), alongside the R packages **ape** and **BioGeoBEARS**.

---

## 📊 Workflow Overview

The entire analysis is structured into three primary stages: a hybrid tree construction, a two-step dating process, and a dual-approach biogeographic analysis.

### 1. Hybrid Phylogenetic Tree Construction
* **Backbone Tree**: A maximum likelihood tree is first built using the `mitogenomes` dataset with **RAxML-NG** to establish a robust phylogenetic backbone.
* **Reconstruction Tree**: The topology of the backbone tree is then fixed and used as a constraint. The more extensive `barcode` dataset is used with **RAxML-NG** to add remaining taxa and re-estimate branch lengths.

### 2. Two-Step Molecular Dating
* **Backbone Dating**: The backbone tree is dated using **treePL**, calibrated with known divergence times between families. This step leverages the higher-quality signal of the mitogenome data.
* **Reconstruction Dating**: Node ages from the dated backbone are extracted and used as multiple secondary calibration points to date the final reconstruction tree using the `ape::chronos` function in R.

### 3. Comprehensive Biogeographic Analysis
* **Ancestral State Reconstruction**: **PASTML** is used to reconstruct ancestral geographic ranges and estimate dispersal routes across the final time-calibrated tree.
* **Time-Stratified Analysis**: **BioGeoBEARS** is employed to test various biogeographic models (e.g., DEC, DIVALIKE). The best-fit model is then used to calculate time-stratified transition matrices, revealing dispersal dynamics across different geological epochs.

---

## ⚙️ System Requirements and Installation

⚠️ **Crucial**: This pipeline is specifically designed for a **Windows** environment using **Windows Subsystem for Linux 2 (WSL2)**. It will not function correctly on other operating systems without modification.

### 1. System Prerequisites
* Windows 10 or 11
* WSL2 with a Linux distribution installed (Ubuntu 22.04 LTS is recommended)

### 2. WSL Dependencies
Open your WSL terminal (e.g., Ubuntu) and run the following commands to install the required command-line tools.

```bash
# Update package lists
sudo apt-get update

# Install RAxML-NG for phylogenetic inference
sudo apt-get install raxml-ng

# Install PASTML for ancestral state reconstruction
sudo apt-get install python3-pip
pip install pastml

# Install treePL for molecular dating (requires compilation)
sudo apt-get install build-essential # Ensure you have build tools
git clone [https://github.com/blackrim/treePL.git](https://github.com/blackrim/treePL.git)
cd treePL/src && make
```
After installation, verify that each tool is accessible from your WSL command line by running commands like `which raxml-ng`.

### 3. R Package Dependencies
Launch R or RStudio and run the following command to install all necessary packages from CRAN and GitHub.

```R
# Install required packages from CRAN
install.packages(c(
  "ape", "stringr", "seqinr", "phytools", "pheatmap", "tidyverse", 
  "ggtree", "ggplot2", "dplyr", "tidyr", "tools", "snow", 
  "MultinomialCI", "reshape2", "gridExtra", "grid", "utils", "base", 
  "phangorn", "plyr", "devtools"
))

# Install BioGeoBEARS from GitHub
devtools::install_github("nmatzke/BioGeoBEARS")
```

---

## 🚀 How to Use

The entire pipeline is controlled and executed from the `tree_biogeography_version1.R` script.

### Step 1: Load Functions
Set your R working directory to the root of this project folder. Then, load all pipeline functions into your R session by sourcing the script.

```R
source("tree_biogeography_version_1.R")
```

### Step 2: Configure Parameters
In your R session, define the variables that will serve as inputs for the main pipeline function. This is where you specify your target taxa, input files, and analysis settings.

```R
# --- 1. Define related groups ---
sister_group_family <- list()
sister_group_family$Staphylinidae<-c("Silphidae_135","Leiodidae_185","Agyrtidae_185","Hydraenidae_185","Ptiliidae_185")
distantly_group<-list()
distantly_group$Staphylinidae<- c("Scarabaeidae", "Tenebrionidae")

# --- 2. Define input file paths (must be relative to the project root) ---
allsequences_path <- "data/test_staphylinidae_barcodes.fasta"
mitogenomes_path <- "data/test_staphylinidae_mitogenomes.fasta"
timeperiods_filepath <- "data/my_timeperiods_20mya.txt"
dispersal_multipliers_filepath <- "data/my_dispersal_multipliers_20mya.txt"

# --- 3. Set computational resources ---
cpu_threads <- 8
```

### Step 3: Run the Main Function
Execute the main pipeline function using the parameters you defined in the previous step. The function will save all output files to the `results/` directory.

```R
# Pass the configuration variables as arguments to the main function
results <- tree_biogeography_pipeline(
  allsequences_path = allsequences_path,
  mitogenomes_path = mitogenomes_path,
  mito_length = 11000,
  sister_group_family = sister_groups,
  distantly_group = distant_groups,
  threads = cpu_threads,
  timeperiods_filepath = timeperiods_filepath,
  dispersal_multipliers_filepath = dispersal_multipliers_filepath,
  construction_model = "GTR+G"
)
```
⚠️ **Note**: The argument names used in the function call (`distantly_group`, `sister_group_family`, etc.) must exactly match the parameter names defined in the `tree_biogeography_pipeline` function.

---

## 📁 Input File Formats

Correct input file formatting is essential for the pipeline to run successfully.

### FASTA Sequence Naming Convention
All sequences in your FASTA files must have headers that follow this strict format, with components separated by underscores (`_`):

`Type_ID_Family_GeographicRegion`

* **Type**: `MMG` for sequences from mitogenomes, or `GMT` for sequences from DNA barcodes.
* **ID**: A unique identifier for the sequence (e.g., a GenBank accession number).
* **Family**: The taxonomic family of the species.
* **GeographicRegion**: The biogeographic region where the species is found.

**Example:**
```fasta
>GMT_JX123456_Staphylinidae_O
ATGCATGCATGC...
>MMG_NC012345_Silphidae_P
ATGCATGCATGC...
```

---

## 📜 License

This project is licensed under the terms of the **MIT License**. See the `LICENSE` file for more details.


