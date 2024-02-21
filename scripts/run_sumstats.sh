#!/bin/bash

# This script creates all sumstats files for the given input files (singleGWAS and cytokine folder)

# Prerequisite: Please download the folder singleGWAS and cytokines from google Drive into the folder "input_folder". 
# Prerequisite: Run the R script significant.R in order to filter for significant GWAS studies. 

significant_singleGWAS_cytokines="../input_files/significant_cytokines_singleGWAS"
sumstats="../output_sumstats/"

# Check if the folder exists
if [ ! -d "$sumstats" ]; then
  # If the folder does not exist, create it
  mkdir -p "$sumstats"
  echo "Folder $sumstats created."
else
  echo "Folder $sumstats already exists."
fi

# Iterate through each file in the unzipped_files and create a sumstats file for each
for file in "$significant_singleGWAS_cytokines"/*; do

    echo "Creating .sumstats file for: $file"

    # Filename
    filename=$(basename "$file" .tsv.gz)

    # Extract the desired part of the filename
    extracted_part=$(echo "$filename" | cut -d'.' -f1)

    # Output filename
    output_filename="../output_sumstats/output_${extracted_part}"

    # Location of snplist
    snplist="../input_files/w_hm3.snplist"

    python2.7 "../ldsc/munge_sumstats.py" --sumstats "$file" \
    --out $output_filename  \
    --snp SNP \
    --N-col sample_size \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --p pvalue \
    --frq EAF \
    --merge-alleles $snplist
done 

