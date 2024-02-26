#!/bin/bash

# This script creates sumstats files for all given input files (in our case 67 GWAS files).
# Sumstats files are necessary for running ldsc later. 

# Prerequisite: Please download the GWAS files into the folder "input_folder".
# Prerequisite: Please download the snplist w_hm3.snplist from ​​https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ into the folder "input_files". 
# Prerequisite: If you have not done it yet, rename the input files that they have the following format "pmidXXXXXX_YY_zzz.tsv". (See README on github for any questions.)


singleGWAS_cytokines_filtered="../input_files/single_GWAS_cytokines_filtered"
sumstats="../output_sumstats/"

# Check if the folder exists
if [ ! -d "$sumstats" ]; then
  # If the folder does not exist, create it
  mkdir -p "$sumstats"
  echo "Folder $sumstats created."
else
  echo "Folder $sumstats already exists."
fi

# Iterate through each file in the  and create a sumstats file for each
for file in "$singleGWAS_cytokines_filtered"/*; do

    echo "Creating .sumstats file for: $file"

    # Filename
    filename=$(basename "$file" .tsv)
    # Extract the part of the filename without .tsv
    extracted_part=$(echo "$filename" | cut -d'.' -f1)

    # Output filename
    output_filename="../output_sumstats/output_${extracted_part}"
    echo "Sumstats output: $output_filename"
    # Location of snplist
    snplist="../input_files/w_hm3.snplist"

    # Run munge_sumstats.py 
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

