# Archive script! 
# This file is for checking the significance of the given 61 GWAS studies in order to only analyze the significant GWAS studies (LDSC and HDL analysing).
# Important: The tools (ldsc and hdl) need the whole GWAS data and no data filtered for significance. 
# Thus, this script is an archived script!


# Set the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Directory of the GWAS files 
folder_singleG_cyt <- "../input_files/single_GWAS_cytokines"
filelist_singleG_cyt <- list.files(folder_singleG, pattern = "\\.tsv\\.gz$", full.names = TRUE)

# Location of the folder where the included significant studies should be stored
included_folderpath <- "../input_files/significant_single_GWAS_cytokines"

# Check if the folder exists
if (!dir.exists(included_folderpath)) {
  # Create it if not existent
  dir.create(included_folderpath, recursive = TRUE)
  print(paste("Folder", included_folderpath, "created."))
} else {
  print(paste("Folder", included_folderpath, "already exists."))
}

# List with files that stores the included (significant) studies
included_files <- c()

# Go through every file
for (file in filelist_singleG_cyt) {
  data <- read.delim(file)
  
  # Check for NA values in "pvalue" column
  amount_na_rows <- sum(is.na(data$pvalue))
  if(amount_na_rows == nrow(data)){
    next # Next if all p values are equal to NA -> Exclude the study 
  }
  
  # Check if p values are < 5*10e-8. 
  small_pvalues <- sum(data$pvalue < 5*10e-8)
  if(small_pvalues == 0){
    next # If all p values are >= 5*10e-8 -> Exclude the study.
  } 
  
  else {
    subset_data <- data[data$pvalue < 5*10e-8, ] # Create a data frame that only contains significant p values. 
    
    # Check if there is only one SNP that fulfills the condition -> Exclude the study. 
    if(nrow(subset_data) == 1){
      next
    }
    
    # Go through every SNP that meets the condition: pvalue < 5*10e-8 
    for (i in 1:(nrow(subset_data) - 1)) {
      # Get the current SNP
      current_snp <- subset_data[i, ]
      
      if(is.na(current_snp$SNP)){
        next # Skip the SNP if the SNP id is NA 
      }
      
      # Flag to check if a pair of SNPs meeting the condition is found
      pair_found <- FALSE
      
      # Iterate through all the other SNPs
      for (j in (i + 1):nrow(subset_data)) {
        # Get the next SNP (don't compare a SNP to itself -> because this is either way not significant)
        next_snp <- subset_data[j, ]
        
        if(is.na(next_snp$SNP)){
          next # Skip the SNP if the SNP id is NA 
        }
        
        
        if (!is.na(current_snp$chr) && !is.na(next_snp$chr)) {
          if (current_snp$chr != next_snp$chr) {
            pair_found <- TRUE
          } else { # If the SNPs are on the same chromosome: Check if the bp are > 1e6 away -> Include the study
            if (!is.na(current_snp$bp_38) && !is.na(next_snp$bp_38)) {
              if (abs(next_snp$bp_38 - current_snp$bp_38) > 1e6) {
                pair_found <- TRUE
              }
            }
          }
        } # end of if condition
        
        
        # If a pair is found, exit the loop -> Include the study 
        if (pair_found) {
          included_files <- c(included_files, file) # add the study to the list for analyizing it further
          break
        }
        
      }
      
      # If a pair is found, exit the outer loop -> Include the study 
      if (pair_found) {
        break
      }
        
        
    } # end of outer for loop
      
      
  } # end of else condition 
    
}

# Loop through each file path in the list and copy it to the destination folder
for (file_path in included_files) {
  file.copy(file_path, included_folderpath)
}

