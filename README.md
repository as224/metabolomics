# Metabolomics

The main goal of this project is exploring genetic correlations and molecular overlaps in inflammatory and autoimmune diseases to prioritize molecular markers and pathways. 


## LDSC 
The LD Score Regression method, designed to estimate heritability and genetic correlation between different GWAS summary statistics, is implemented as a tool. The main repository for it can be found here: https://github.com/bulik/ldsc.
In our analysis, we applied their algorithm to 67 GWAS studies. The list of their ID's can be found at /input_files/GWAS_ids.tsv. 

### Setup Instructions for LDSC

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

Clone the proper LDSC GitHub repository ([https://github.com/bulik/ldsc.git](https://github.com/bulik/ldsc.git)) into the metabolomics repository with the following structure: `metabolomics/ldsc`.

Following these steps will organize the data effectively for using this repository properly.

Below one can familiarize with our workflow for utilizing LDSC and HDL. It is crucial to follow these steps in order to reproduce our results.

## Workflow of LDSC and HDL
![Alt text](https://github.com/as224/metabolomics/blob/main/workflow_ldsc_hdl.png "Workflow LDSC and HDL")


