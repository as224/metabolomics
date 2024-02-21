#!/bin/bash

# Path to the folder containing the sumstats files
sumstats="../output_sumstats/" # The folder should already exist (look at run_sumstats.sh if not)
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
    file_list="${file1},${file1},"
    amount_correlation_files=1

    # Iterate through the files again to check for unique middle parts
    for ((j = 0; j < ${#files[@]}; j++)); do
        file2="${files[j]}" # Second file in the list 

        # Extract the filename without extension for file2
        filename2=$(basename "$file2" | cut -d'.' -f1)

        # If file1 and file2 are the same, skip
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
    echo "The amount of files that we've correlated: $amount_correlation_files"

    # Output filename
    output_filename="../output_ld_regression/${filename1}_ldReg"

    # Run ldsc.py for file1 and file2. 
    # Prerequisites 1: Please download the eur_w_ld_chr folder from https://zenodo.org/records/8182036 and put it into the folder "input_files"    output_filename="../output_ld_regression/${filename1}_ldReg"

    python2.7 "../ldsc/ldsc.py" --rg $file_list \
     --ref-ld-chr "../input_files/eur_w_ld_chr/" \
     --w-ld-chr "../input_files/eur_w_ld_chr/" \
     --out $output_filename
    
done



