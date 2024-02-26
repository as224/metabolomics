#!/bin/bash

# This script aims to run the script adjustingEverything.py (check column length, include sample size and cut SNPs that have EAF=. or =NA)

# python_script should be in the same folder as this script 
python_script="adjustingEverything.py"

input_file_path="../input_files/single_GWAS_cytokines"
output_file_path="../input_files/single_GWAS_cytokines_complete"
output_file_path_complete="../input_files/single_GWAS_cytokines_filtered"

# Check if folder exists, otherwise create it 
if [ ! -d "$output_file_path" ]; then
  # If the folder does not exist, create it
  mkdir -p "$output_file_path"
  echo "Folder $output_file_path created."
else
  echo "Folder $output_file_path already exists."
fi

# Check if folder exists, otherwise create it 
if [ ! -d "$output_file_path_filtered" ]; then
  # If the folder does not exist, create it
  mkdir -p "$output_file_path_filtered"
  echo "Folder $output_file_path_filtered created."
else
  echo "Folder $output_file_path_filtered# already exists."
fi

for file in "$input_file_path"/*; do    
  filename=$(basename -- "$file")
  
  output_filename="${output_file_path}/${filename}"
  output_filename_filtered="${output_file_path_filtered}/${filename}"
  
  echo "Output filename (same amount of columns per SNP):${output_filename}"
  echo "Output filename (filtered for eaf and sample size included where missing):${output_filename_filtered}"

  # run python script  
  python3.9 "$python_script" "$file" "$output_filename" "$output_filename_filtered"
done


