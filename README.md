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

### Requirements and Workflow:
The chosen GWAS study that is harmonized, should contain the follwing values for all investigated Single Nucleotide Polymorphisms (SNPs):
- SNP ID
- Chromosome on which the SNP is located
- Genomic Position or base pair coordinate, preferably on assembly GRCH38 (otherwise, transformation is later possible)
- Effect allele
- Non-effect-allele
- Effect Allele Frequency (EAF), representing the frequency of the effect allele in the population
- Effect direction or beta effect
- Odds Ratio (OR), representing the association between the SNP and the trait
- Standard Error (SE)
- P-value
- Investigated trait or phenotype
- Sample size of the study


Starting, run Harmonize.sh. Here, adapt the script to the available columns of the GWAS file. The resulting TSV contains the necessary columns in the needed format. 


If the data is from an earlier assembly, the tool "LiftOver" can be used for transformation. For that, add the --bed parameter, and the script generates the necessary input file for LiftOver. Following, the LiftOver results are merged with the origin file to obtain the new positions: 
```bash
awk -F"\t" 'FNR==NR{a[$4]=$2; next}{print $1"\t"$2"\t"a[$1]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' liftOverResult.bed resultHarmonize.tsv > result.final
```

Furthermore, the script can filter p-values and SNPs with unavailable EAFs if needed. 

Following, run 2_merger_colnames.pl with the resulting file. This script switches the effect and non-effect alleles if the EAF is over 0.5 and adapts the EAF, beta effect, and odds ratio (OR) accordingly. If dbSNP is available, allele frequencies and/or the Minor Allele Frequency (MAF) can be added. 

### Scripts:

#### 1. Harmonize.sh
```bash
bash harmonize.sh -i gwas -t trait [--pvalue] [--bed] [--eaf]
```


Parameters:
- -i: The GWAS study that is harmonized (in gzip format)
- -t: the studied trait
- --pvalue: include, if SNPs with a p-value >= 0.05 should be excluded
- --bed: include, if a .bed file should be created
- --eaf: include, if invalid EAFs (NA or “.”) should be filtered


#### 2. 2_merger_colnames.pl
```perl
perl 2_merger_colnames.pl resulting_SH_file dbsnp output_file
```


Parameters:
- resulting_SH_file: output file of bring_to_gwascatalog_format.sh
- dbsnp: theoretically, here, dbsnp files can be used to add allele frequency from dbSNP to the study file and/or fix MAF in case it was provided.
- output_file: output directory



#### 3. other scripts
- `scripts_harmonize/bring_to_gwascatalog_format.sh`: Second origin script for the second steps of Harmonization while creating the .bed files.
- `scripts_harmonize/preprocess.sh`: First Origin script for the first step of the Harmonization
- `scripts_harmonize/filter_pval.py`: alternative script for p-value filtering
- `scripts_harmonize/filter_EAF.py`: used during Harmonization if --pvalue is added, filters SNPs that have an invalid EAF



## LDSC <a name="LDSC-overall"></a>
The LD Score Regression method, designed to estimate heritability and genetic correlation between different GWAS summary statistics, is implemented as a tool. The main repository for it can be found here: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc). There one can also find the necessary packages and versions in order to run ldsc. The original paper is linked here [LD Score regression distinguishes confounding from polygenicity in genome-wide association studies](https://doi.org/10.1038/ng.3211).
In our analysis, we applied their algorithm to 67 GWAS studies. The list of their ID's can be found at `input_files/GWAS_ids.tsv`. 

### Setup Instructions for LDSC <a name="setup-instructions-ldsc"></a>

To maintain the proper folder structure in your cloned GitHub repository, please follow these steps:

#### 1. Create a Folder for Input Files

Please create a new folder named "single_GWAS_cytokines" where you will store all the input files.

#### 2. Download GWAS Studies

Download all the desired GWAS studies and save them as input files in the `metabolomics/input_files` folder. The filenames should adhere to the following structure:

- **File Naming Convention**: pmid00000000_XX_eur.tsv
  - `pmid00000000`: Replace zeros with the specific PubMed ID number of the study.
  - `XX`: Specification of the disease studied.
  - `eur` (or `eas`, etc.): Indicates the continent of the study's population.

For example:
- `pmid21833088_MS_eur.tsv`
- `pmid29632382_T2d-BMIadj_eur.tsv`

#### 3. Download Additional Files

In the `metabolomics/input_files` folder, also download the following files from the provided sources:

- **w_hm3.snplist:**
  - Download from: [w_hm3.snplist](https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist)
  - Save in `metabolomics/input_files`

- **eur_w_ld_chr Folder:**
  - Download from: [eur_w_ld_chr Folder](https://zenodo.org/records/8182036)
  - Save the entire folder here, that the structure looks like that: `metabolomics/input_files/eur_w_ld_chr`

#### 4. Clone the Proper LDSC GitHub Repository

Clone the proper LDSC GitHub repository ([https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)) into the metabolomics repository with the following structure: `metabolomics/ldsc`.

Following these steps will organize the data effectively for using this repository properly.


## HDL <a name="HDL-overall"></a>
As an extension to LDSC, the High-Definition Likelihood (HDL) method was developed to enhance precision in estimating genetic correlations. The tool generates summary statistics for input GWAS files and calculates the genetic correlation between GWAS pairs. 
For a more comprehensive understanding of the tool, detailed information can be found on the associated GitHub repository: [https://github.com/zhenin/HDL](https://github.com/zhenin/HDL). The original paper is linked here [High-definition likelihood inference of genetic correlations across human complex traits](https://doi.org/10.1038/s41588-020-0653-y).
In our analysis, we also employed the HDL method to the 67 GWAS studies mentioned above. `input_files/GWAS_ids.tsv`

### Setup Instructions for HDL <a name="setup-instructions-hdl"></a>

To maintain the proper folder structure in your cloned GitHub repository, please follow these steps:

#### 1. Create a Folder for Input Files and download GWAS Studies

Follow the instructions in steps 1 and 2 at [Setup Instructions for LDSC](#setup-instructions-ldsc).

#### 2. Download Additional Files

In the `metabolomics/input_files` folder, also download the following folder from the provided source:

- **UKB_array_SVD_eigen90_extraction Folder:**
  - Download from: [UKB_array_SVD_eigen90_extraction](https://www.dropbox.com/s/fuvpwsf6r8tjd6c/UKB_array_SVD_eigen90_extraction.tar.gz?dl=0)
  - Save the entire folder "UKB_array_SVD_eigen90_extraction" in the `metabolomics/input_files` folder.
  - You can use another refence panel here, if needed. Just adjust the path in run_wrangling.R and run_hdl.R.

#### 3. Clone the Proper HDL GitHub Repository

Clone the proper HDL GitHub repository ([https://github.com/zhenin/HDL](https://github.com/zhenin/HDL)) into the metabolomics repository with the following structure: `metabolomics/HDL`.

Following these steps will organize the data effectively for using this repository properly.

In order to understand our analysis, please follow the [Workflow of LDSC and HDL](#workflow-ldsc-hdl).


## Workflow of LDSC and HDL <a name="workflow-ldsc-hdl"></a>
The following figure describes our workflow for utilizing the LDSC and HDL approach. It is crucial to follow these steps in order to reproduce our results. 

![Alt text](https://github.com/as224/metabolomics/blob/main/workflow_ldsc_hdl.png "Workflow LDSC and HDL")


