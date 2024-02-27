# Metabolomics

The main goal of this project is exploring genetic correlations and molecular overlaps in inflammatory and autoimmune diseases to prioritize molecular markers and pathways. 

## Table of Contents
1. [Harmonization of a GWAS Study](#harmonize)
2. [LDSC](#LDSC-overall)
3. [Setup Instructions for LDSC](#setup-instructions-ldsc)
4. [HDL](#HDL-overall)
5. [Setup Instructions for HDL](#setup-instructions-hdl)
6. [Workflow of LDSC and HDL](#workflow-ldsc-hdl)

## Harmonization <a name="harmonize"></a>
The main objective here is to ensure a consistent format for the GWAS studies utilized in the subsequent analytical procedures. Consequently, the availability of the necessary values is guaranteed.

### Requirements:

### 1. Harmonize.sh
Command: bash harmonize.sh -i gwas -t trait [--pvalue] [--bed] [--eaf] 

Parameters: 
- -i: The GWAS study that is harmonized (in gzip format)
- -t: the studied trait
- --pvalue: include, if SNPs with a p-value >= 0.05 should be excluded
- --bed: include, if a .bed file should be created
- --eaf: include, if invalid EAFs should be filtered

### 2. 2_merger_colnames.pl
Command: perl 2_merger_colnames.pl  resulting_SH_file  dbsnp output_file

Parameters: 
- resulting_SH_file: output file of bring_to_gwascatalog_format.sh
-	dbsnp: theoretically, here, dbsnp files can be used to add allele frequency from dbSNP to the study file and/or fix MAF in case it was provided.
-	output_file: output directory


### 3. other scripts
- bring_to_gwascatalog_format.sh: 
- preprocess.sh:
- filter_pval.py: alternative script for p-value filtering
- filter_EAF.py: used during Harmonization if --pvalue is added



## LDSC <a name="LDSC-overall"></a>
The LD Score Regression method, designed to estimate heritability and genetic correlation between different GWAS summary statistics, is implemented as a tool. The main repository for it can be found here: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc).
In our analysis, we applied their algorithm to 67 GWAS studies. The list of their ID's can be found at /input_files/GWAS_ids.tsv. 

### Setup Instructions for LDSC <a name="setup-instructions-ldsc"></a>

To maintain the proper folder structure in your cloned GitHub repository, please follow these steps:

#### 1. Create a Folder for Input Files

Please create a new folder named "single_GWAS_cytokines" where you will store all the input files.

#### 2. Download GWAS Studies

Download all the desired GWAS studies and save them as input files in the "input_files" folder. The filenames should adhere to the following structure:

- **File Naming Convention**: pmid00000000_XX_eur.tsv
  - `pmid00000000`: Replace zeros with the specific PubMed ID number of the study.
  - `XX`: Specification of the disease studied.
  - `eur` (or `eas`, etc.): Indicates the continent of the study's population.

For example:
- `pmid21833088_MS_eur.tsv`
- `pmid29632382_T2d-BMIadj_eur.tsv`

#### 3. Download Additional Files

In the "input_files" folder, also download the following files from the provided sources:

- **w_hm3.snplist:**
  - Download from: [w_hm3.snplist](https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist)
  - Save in the "input_files" folder.

- **eur_w_ld_chr Folder:**
  - Download from: [eur_w_ld_chr Folder](https://zenodo.org/records/8182036)
  - Save the entire folder "eur_w_ld_chr" in the "input_files" folder.

#### 4. Clone the Proper LDSC GitHub Repository

Clone the proper LDSC GitHub repository ([https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)) into the metabolomics repository with the following structure: `metabolomics/ldsc`.

Following these steps will organize the data effectively for using this repository properly.


## HDL <a name="HDL-overall"></a>
As an extension to LDSC, the High-Definition Likelihood (HDL) method was developed to enhance precision in estimating genetic correlations. The tool generates summary statistics for input GWAS files and calculates the genetic correlation between GWAS pairs. 
For a more comprehensive understanding of the tool, detailed information can be found on the associated GitHub repository: [https://github.com/zhenin/HDL](https://github.com/zhenin/HDL).
In our analysis, we also employed the HDL method to the 67 GWAS studies mentioned above. (/input_files/GWAS_ids.tsv)

### Setup Instructions for HDL <a name="setup-instructions-hdl"></a>

To maintain the proper folder structure in your cloned GitHub repository, please follow these steps:

#### 1. Create a Folder for Input Files and download GWAS Studies

Follow the instructions in steps 1 and 2 at [Setup Instructions for LDSC](#setup-instructions-ldsc).

#### 2. Download Additional Files

In the "input_files" folder, also download the following folder from the provided source:

- **UKB_array_SVD_eigen90_extraction Folder:**
  - Download from: [UKB_array_SVD_eigen90_extraction](https://www.dropbox.com/s/fuvpwsf6r8tjd6c/UKB_array_SVD_eigen90_extraction.tar.gz?dl=0)
  - Save the entire folder "UKB_array_SVD_eigen90_extraction" in the "input_files" folder.
  - You can use another refence panel here, if needed. Just adjust the path in run_wrangling.R and run_hdl.R.

#### 3. Clone the Proper HDL GitHub Repository

Clone the proper HDL GitHub repository ([https://github.com/zhenin/HDL](https://github.com/zhenin/HDL)) into the metabolomics repository with the following structure: `metabolomics/HDL`.

Following these steps will organize the data effectively for using this repository properly.

In order to understand our analysis, please follow the [Workflow of LDSC and HDL](#workflow-ldsc-hdl).


## Workflow of LDSC and HDL <a name="workflow-ldsc-hdl"></a>
The following figure describes our workflow for utilizing the LDSC and HDL approach. It is crucial to follow these steps in order to reproduce our results. 

![Alt text](https://github.com/as224/metabolomics/blob/main/workflow_ldsc_hdl.png "Workflow LDSC and HDL")


