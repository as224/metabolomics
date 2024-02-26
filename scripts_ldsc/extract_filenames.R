# This script is for extracting the filenames in order to get a list of all the GWAS we analyzed. 

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

folder_path <- "../input_files"

# Set input path and create list of included files 
folder_singleG <- "../input_files/single_GWAS_cytokines_filtered"
filelist_singleG_cyt <- list.files(folder_singleG, pattern = "\\.tsv$", full.names = TRUE)

# Pattern to match
pattern <- "pmid(\\d+)_(.+?)\\.tsv"
result_list <- vector("character", length(filelist_singleG_cyt))

# Extract matches
matches <- gregexpr(pattern, filelist_singleG_cyt, perl = TRUE)
matches_list <- regmatches(filelist_singleG_cyt, matches)

# Loop through matches
for (i in seq_along(matches_list)) {
  print(matches_list[[i]])
  sub_matches <- unlist(strsplit(matches_list[[i]], "_"))
  continent_tsv <- sub_matches[3]
  continent <- unlist(strsplit(continent_tsv, "\\."))
  
  
  pmid <- sub_matches[1]
  condition <- paste(sub_matches[2], "_", continent, sep = "")
  result_list[i] <- paste(pmid, "(", condition, ")", sep = "")
  
  
}

result_list_clean <- gsub("pmid|\\.tsv\\.gz", "", result_list)
result_str <- paste(result_list_clean, collapse = ", ")

result_df <- as.data.frame(result_str)

# Print the result
write.table(result_df, file = paste(folder_path, "/GWAS_ids.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
 