#!/bin/bash

# Path to the folder containing the sumstats files
# Prerequisites: Please download the eur_w_ld_chr folder from https://zenodo.org/records/8182036 and put it into the folder "input_files"   

sumstats="../output_sumstats/" 
output_ldsc="../output_ld_regression/"

if [ ! -d "$output_ldsc" ]; then
  # If the folder does not exist, create it
  mkdir -p "$output_ldsc"
  echo "Folder $output_ldsc created."
else
  echo "Folder $output_ldsc already exists."
fi

# Take only the sumstats files that end with ".gz". This is necessary for running ldsc.py. 
files=("$sumstats"*.gz)

for ((i = 0; i < ${#files[@]}; i++)); do
    file1="${files[i]}" # First file in the list (as input for ldsc.py)
    filename1=$(basename "$file1" | cut -d'.' -f1) # Extract the filename without extension
    
    # Create a string that stores all the filenames 
    file_list="${file1},${file1}," # The first two elements of the list are the same file in order to calculate the correlation of GWAS1-GWAS1. This genetic correlation value should be 1. 
    amount_correlation_files=1

    # Iterate through the files again to add them to the correlation list 
    for ((j = 0; j < ${#files[@]}; j++)); do
        file2="${files[j]}" 

        # Extract the filename without extension for file2
        filename2=$(basename "$file2" | cut -d'.' -f1)

        # If file1 and file2 are the same, skip (as we handeled this above)
        if [ "$file1" == "$file2" ]; then
            continue
        fi

        # Add file to the list 
        file_list="${file_list}${file2},"
        ((amount_correlation_files++))

    done
    #Remove the trailing comma
    file_list=${file_list%,}

    echo "File list for $filename1: $file_list"
    echo "The amount of correlated files: $amount_correlation_files"

    # Output filename
    output_filename="../output_ld_regression/${filename1}_ldReg"

    # Run ldsc.py for file1 and file2. 
    python2.7 "../ldsc/ldsc.py" --rg $file_list \
     --ref-ld-chr "../input_files/eur_w_ld_chr/" \
     --w-ld-chr "../input_files/eur_w_ld_chr/" \
     --out $output_filename
    
done



