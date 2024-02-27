# To rename the single_GWAS full files
# Example:
# pmid21833088_MS_eur_full.merge.gz
# to pmid21833088_MS_eur.tsv

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

input_directory <- "../input_files/single_GWAS_full"

# Files to rename
files <- list.files(input_directory)

# Rename files
for (file in files) {
  if (grepl("\\.gz$", file)) {
    # Old file path (files to rename)
    old_file_path <- paste(input_directory, "/", file, sep="")
    
    # Get new file name
    new_name <- paste(sub("^((?:[^_]*_){2}[^_]*).*", "\\1", file), ".tsv", sep = "") 
    new_file_path <- paste(input_directory, "/", new_name, sep="")
    
    # Rename old file 
    file.rename(old_file_path, new_file_path)
  }
}

