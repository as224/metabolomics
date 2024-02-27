# To filter out the hdl sumstats that reach min. 95% overlap with reference panel 
# HDL prints a warning: "Warning: More than 1% SNPs in reference panel are missed in GWAS"

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#install.packages("readtext")
library("readtext")

input_directory <- "../output_wrangling" # path to wrangling output

# Minimum overlap we want between gwas and reference panel
min_overlap <- 95 # Adjust as needed

files <- list.files(path = input_directory, pattern = "\\.txt$", full.names = TRUE)

# Loop over directory and remove files 
for (txt_file in files) {
  # Read file
  log <- readLines(txt_file)
  
  # Loop over lines and find the percentage line
  for (line in log) {
    if (grepl("SNPs in reference panel are available in GWAS.", line)) {
      # Extract the percentage using regex
      percentage <- regmatches(line, regexpr("\\((.*?)%\\)", line))[[1]]
      # Remove the parentheses and percent sign
      percentage <- gsub("\\(|\\)|%", "", percentage)
      
      # If percentage is below 95% remove the GWAS (txt and rds file)
      if (percentage < min_overlap) {
        # Get rds file path
        rds_file <- sub("\\.txt$", ".hdl.rds", txt_file)
        # Remove files
        file.remove(txt_file)
        file.remove(rds_file)
      }
    } 
  }
}