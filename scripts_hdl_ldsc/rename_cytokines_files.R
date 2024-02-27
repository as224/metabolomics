# To rename the cytokines full files
# Example: 
# bNGF_27989323-GCST004421-EFO_0008035.h.tsv.merge.full
# to pmid27989323_bNGF_eur_full.tsv

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

input_directory <- "../input_files/cytokines_full"

# Files to rename
files <- list.files(input_directory)

# Rename files
for (file in files) {
  # Old file path (files to rename)
  old_file_path <- paste(input_directory, "/", file, sep="")
  print(old_file_path)
  # Split file name into parts 
  parts <- unlist(strsplit(file, "[-_]"))
  
  # Extract parts for new file name
  pmid <- parts[2] # numerical pmid 
  cytokine <- parts[1] # trait name
  region <- "eur" # region
  
  # Get new file name
  new_name <- paste("pmid", pmid, "_", cytokine, "_", region, ".tsv", sep="")
  new_file_path <- paste(input_directory, "/", new_name, sep="")
  
  # Rename old file 
  file.rename(old_file_path, new_file_path)
}