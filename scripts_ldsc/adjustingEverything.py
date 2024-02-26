# This script is for adjusting the given GWAS files in order to be able to use ldsc and hdl. 

# Import packages
import csv
import sys
import pandas as pd 
import os 

# Set working directory 
working_dir = os.path.dirname(os.path.abspath(__file__))

# The aim of this function is to only keep the SNPs of a GWAS where the amount of columns are equal the amount of headers of the file.
def process_colLength(input, output_file_path):   
    print("PROCESSING: ", input)
    
    # Open tsv file to process it
    with open(input, 'r', newline='', encoding='utf-8') as file:
        count = 0 
        header = file.readline().strip().split('\t') # header row 
        print("Amount of columns: ", len(header)) # length of header row 

        # Open new file 
        with open(output_file_path, 'w', newline='', encoding='utf-8') as output:
            tsv_writer = csv.writer(output, delimiter='\t')
            tsv_writer.writerow(header)

            for line in file:
                    columns = line.strip().split('\t')
                    
                    # Check if amount of entries for a column is equal to the amount of headers
                    if len(columns) == len(header):
                        output.write(line) # If equal -> keep the line
                    else:
                        count += 1 # If unequal -> cut the line out and count the count + 1 in order to see how many lines have been cut out 

    print("Amount of rows with column length bigger than", len(header), ": ", count)
    print("Saved: ", output_file_path)

# The aim of this function is to only keep the SNPs that have a real EAF value (no . or NA) and replace the sample size with the real sample size if it is NA 
def process_eafsAndsampleSize(input_file_path, output_file_path):
    # Get the sample sizes to the according GWAS studies
    matching_dataframe_path = "../input_files/sample_size_cytokines.tsv"
    sample_size_df = pd.read_csv(matching_dataframe_path, sep="\t")
   
    print("PROCESSING pre-filtered output file: ", input_file_path)
    
    # Get filename
    input_filename = os.path.basename(input_file_path)
    # Read the input file 
    df_input_filtered = pd.read_csv(input_file_path, sep="\t")
    
    # Check if the input file starts with pmid27989323 -> in these files the sample sizes are missing
    if input_filename.startswith("pmid27989323"):
        result = sample_size_df['id'] == input_file_path # Check which filename fits 
        sample_size = sample_size_df.loc[result, 'sample_size'].values[0] if any(result) else None # Get sample size 
        df_input_filtered['sample_size'] = sample_size # Add the sample size 
        
    # Drop all rows that contain NA or . for EAF 
    rows_before_EAFdrop = len(df_input_filtered)
    df_input_filtered['EAF'] = df_input_filtered['EAF'].replace("NA", pd.NA)
    df_input_filtered['EAF'] = df_input_filtered['EAF'].replace(".", pd.NA)
    df_input_filtered = df_input_filtered.dropna(subset=['EAF'])
    rows_after_EAFdrop = len(df_input_filtered)
    # Get the amount of dropped rows (because of EAF=. or =NA)
    rows_dropped = rows_before_EAFdrop - rows_after_EAFdrop
    
    # Pvalue and sample size column as numeric
    df_input_filtered['pvalue'] = pd.to_numeric(df_input_filtered['pvalue'], errors='coerce')
    df_input_filtered['sample_size'] = pd.to_numeric(df_input_filtered['sample_size'], errors='coerce')
    
    # Save df_input_filtered as csv
    df_input_filtered.to_csv(output_file_path, sep='\t', index=False)
    
    print("Total amount of dropped rows: ", rows_dropped)
    print("Ergebnisse in", output_file_path, "gespeichert.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Parameters are missing!")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    output_file_complete = sys.argv[3]
    process_colLength(input_file, output_file)
    process_eafsAndsampleSize(output_file, output_file_complete)

    
   