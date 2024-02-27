# To extract the sample size from old cytokines files (significance filtered) because they are missing in the full files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sample_info_path <- "../input_files_significant/cytokines" # files containing the sample size
output_tsv <- "../input_files/sample_size_cytokines.tsv" # output tsv (id, sample_size)

# Find the files to open
files <- list.files(sample_info_path)

# Create an empty data frame to store all sample sizes
column_names <- c("id", "sample_size") # col names
sample_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(sample_df) <- column_names

# Rename files
for (file in files) {
  # Read file
  file_path <- paste(sample_info_path, "/", file, sep = "")
  data <- head(read.table(gzfile(file_path), header = TRUE, sep = "\t"))
  selected_columns <- data[1:1,] # Read first two rows
  
  # Create empty df for file
  cols <- c("id", "sample_size")
  df <- data.frame(matrix(ncol = length(cols), nrow = 1))
  colnames(df) <- cols
  
  # Fill cells 
  id <- paste("../input_files/single_GWAS_cytokines_complete/", sub("\\.tsv.*", ".tsv", file), sep = "")
  df$id[1] <- id
  df$sample_size[1] <- selected_columns$sample_size
  
  # Add file and sample_size to sample_df
  sample_df <- rbind(sample_df, df)
}

# Write final df to tsv
write.table(sample_df, file = output_tsv, sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)

