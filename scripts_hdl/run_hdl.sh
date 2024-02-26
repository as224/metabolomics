#!/bin/bash

# Path to the folder containing the sumstats files
sumstats="../output_wrangling/" # The folder should already exist (look at run_wrangling.sh if not)
output_hdl="../output_hdl/"

if [ ! -d "$output_hdl" ]; then
  # If the folder does not exist, create it
  mkdir -p "$output_hdl"
  echo "Folder $output_hdl created."
else
  echo "Folder $output_hdl already exists."
fi

# Array to store files with ".rds" extension
files=("$sumstats"*.rds)

# Loop through the array to create all possible combinations
for ((i = 0; i < ${#files[@]}; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do
        file1="${files[$i]}" # First file in the list (as input for HDL.run.R)
        filename1=$(basename "$file1" | cut -d'.' -f1) # Extract the filename without extension
        
        file2="${files[$j]}"
        filename2=$(basename "$file2" | cut -d'.' -f1) # Extract the filename without extension 

        output_filename="../output_hdl/${filename1}_${filename2}.raw.gwas.Rout"
        
        # Print gwas pair
        echo "Applying HDL on $file1 and $file2"

        # Command
        Rscript "../HDL/HDL.run.R" \
        gwas1.df=$file1 \
        gwas2.df=$file2 \
        LD.path="../reference_panels_hdl/UKB_array_SVD_eigen90_extraction" \
        output.file=$output_filename

    done
done




