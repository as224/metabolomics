#!/bin/bash

# This script creates all hdl result files (hdl(gwas1,gwas2)) for the given input sumstats files (output_wrangling)

# Prerequisite: Please download the folder singleGWAS and cytokines from google Drive into the folder "input_files". 
# Prerequisite: Run run_adjustingEverything.sh before to check and adjust the tsv files
# Prerequisite: Run the run_wrangling.sh script to generate the sumstats files

# Path to the folder containing the sumstats files
sumstats="../output_wrangling/" # The folder should already exist (look at run_wrangling.sh if not)
output_hdl="../output_hdl/" # output path

if [ ! -d "$output_hdl" ]; then
  # If the folder does not exist, create it
  mkdir -p "$output_hdl"
  echo "Folder $output_hdl created."
else
  echo "Folder $output_hdl already exists."
fi

# Array to store files with ".rds" extension
files=("$sumstats"*.rds)

# Loop through the array to create file pairs where the files are different and avoid counting the reverse order
for ((i = 0; i < ${#files[@]}; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do  # Start the inner loop at i + 1 to avoid counting the reverse order
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
        LD.path="../input_files/UKB_array_SVD_eigen90_extraction" \
        output.file=$output_filename \
        numCores=4

    done
done




