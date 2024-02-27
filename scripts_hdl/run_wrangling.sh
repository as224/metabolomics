#!/bin/bash

# This script creates all sumstats files for the given input files (singleGWAS and cytokine folder)

# Prerequisite: Please download the folder singleGWAS and cytokines from google Drive into the folder "input_files". 
# Prerequisite: Run run_adjustingEverything.sh before to check and adjust the tsv files

significant_singleGWAS_cytokines="../input_files/single_GWAS_cytokines_filtered" # GWAS input files
sumstats="../output_wrangling/" # output path

# Check if the sumstats folder exists
if [ ! -d "$sumstats" ]; then
  # If the folder does not exist, create it
  mkdir -p "$sumstats"
  echo "Folder $sumstats created."
else
  echo "Folder $sumstats already exists."
fi

# Iterate through each file in the single_GWAS_cytokines_filtered and create a sumstats file for each
for file in "$significant_singleGWAS_cytokines"/*; do

    echo "Creating .sumstats file for: $file"

    # Filename
    filename=$(basename "$file" .tsv.gz)

    # Extract the desired part of the filename
    extracted_part=$(echo "$filename" | cut -d'.' -f1)

    # Output filename
    output_filename="../output_wrangling/output_${extracted_part}"

    # Command
    Rscript "../HDL/HDL.data.wrangling.R" \
    gwas.file="$file" \
    LD.path=../input_files/UKB_array_SVD_eigen90_extraction \
    SNP=SNP A1=effect_allele A2=non_effect_allele N=sample_size b=beta_effect se=SE \
    output.file="$output_filename" \
    log.file="$output_filename"
done 
