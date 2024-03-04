# Archive script! 
# Script to calculate the SNP overlap between the two GWAS sets cytokines and single_GWAS

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function to calculate the percentage overlap
overlap <- function(gwasA, gwasB) {
  overlap <- length(intersect(gwasA$SNP, gwasB$SNP))
  print(overlap)
  perc_overlap <- (overlap/length(gwasA$SNP))*100
  return(perc_overlap)
}

# Files
file_list_cytokines <- list.files(path = "../input_files/cytokines", pattern = "\\.tsv.gz$", full.names = TRUE)
file_list_GWAS <- list.files(path = "../input_files/single_GWAS", pattern = "\\.tsv.gz$", full.names = TRUE)
file_list_all <- c(file_list_cytokines, file_list_GWAS)

# Generate all possible file pairs
file_pairs <- combn(file_list_all, 2)

# Iterate over pairs
overlaps <- list()
for (i in seq(ncol(file_pairs))) {
  file1 <- file_pairs[1, i]
  file2 <- file_pairs[2, i]
  
  # Load files
  gwasA <- read.table(gzfile(file1), header = TRUE, sep = "\t") 
  gwasB <- read.table(gzfile(file2), header = TRUE, sep = "\t")  
  
  # Apply function
  result <- overlap(gwasA, gwasB)
  
  # Print overlap
  print(result)
  
  # Save results
  overlaps[[i]] <- result
}
