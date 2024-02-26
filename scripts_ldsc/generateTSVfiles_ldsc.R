# The goal of this script is to generate tsv files including all the important information of the ldsc output. 
# Also, the aim is to generate correlation matrices with real, corrected and significant correlation values.

# Load required libraries
# install.packages("readtext")
library("readtext")
library(tidyr)
library(dplyr)
library(tibble)

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set the path to where the output log files of ldsc are stored
output_ld_regression_path <- "../output_ld_regression"

# Function that generates tsv from log and adds cols (p1_id, p2_id, id and method)
generate_df <- function(input) {
  log <- readLines(input)
  
  # Start and end of table (of the log file)
  start_index <- which(log == "Summary of Genetic Correlation Results")
  end_index <- grep("^Analysis finished at", log)
  
  # Extract only the interessting part of the log file 
  cut_log <- log[(start_index + 1):(end_index - 2)]
  
  # Split the lines into separated cols (using spaces)
  data <- strsplit(cut_log, "\\s+") 
  
  # Store the data in a data frame
  df <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
  colnames(df) <- df[1, ] # Set col names
  df <- df[-1,-1]  # Delete first row (names) and col (empty) 
  row.names(df) <- NULL # Reset indices
  
  # Add ID column
  # Extract p1 and p2
  p1_p2 <- df[, 0:2]
  
  extract_substring <- function(string) {
        sub("^(?:[^_]*_){2}", "", gsub(".sumstats.gz$", "", string))
      }
  # Apply to p1 and p2
  p1_p2 <- apply(p1_p2, 2, extract_substring)
  
  # Add four additional cols to df
  df$p1_id <- p1_p2[,1] # ID of gwas 1
  df$p2_id <- p1_p2[,2] # ID of gwas 2
  df$id <- paste(df$p1_id, df$p2_id, sep = "/") # Matching ID (ID1/ID2)
  df$method <- rep("ldsc", nrow(df)) # Used method
  
  # Replace NA with 0 for genetic correlation 
  df$rg <- as.numeric(df$rg)
  df$rg[is.na(df$rg)] <- 0

  # Return the df
  return(df)
}

# Initialize an empty data frame without specifying the number of rows
df_all <- data.frame(matrix(ncol = 16, nrow=0))
colnames(df_all) <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se", "p1_id", "p2_id", "id", "method")

files <- list.files(path = output_ld_regression_path, pattern = "\\ldReg.log$", full.names = TRUE)
counter_files=0

# Create a tsv file for every output log file 
for (file in files) {
  input_file_path <- file
  output_file_path <- gsub(".log", ".tsv", input_file_path)

  # Call function
  df <- generate_df(input = input_file_path)
  counter_files = counter_files+1

  # Write df to tsv file
  write.table(df, file = output_file_path , sep = "\t", row.names = FALSE, quote=FALSE)
  
  # Append df to df_all
  df_all <- rbind(df_all, df)
}

# Adjust the correlation values -> >1=0 and <-1=0 (the highest correlation (diagnoal values) is one -> thus a higher value than 1 explains something weird and also <-1)
df_all_adjusted <- df_all  
df_all_adjusted$rg[df_all$rg > 1] <- 0
df_all_adjusted$rg[df_all$rg < -1] <- 0
df_all_adjusted$p <- as.numeric(df_all_adjusted$p)
df_all_adjusted$q <- p.adjust(df_all_adjusted$p, method = "fdr")

# Dataframe with corrected and significant p-values
df_all_adjusted_significant <- df_all_adjusted[df_all_adjusted$p < 0.05 & !(is.na(df_all_adjusted$p)), ]

# Create a matrix containing all correlation values (with adjusted correlation values -> >1, <-1 handeled)
df_aggregated <- df_all_adjusted %>%
  group_by(p1_id, p2_id) %>%
  summarise(rg = first(rg)) %>%
  ungroup()

matrix_df <- spread(df_aggregated, key = p1_id, value = rg)

row_names <- matrix_df$p2_id
matrix_df$p2_id <- NULL
rownames(matrix_df) <- row_names

matrix_df <- as.matrix(matrix_df)

# Create a matrix containing all real (not adjusted) correlation values
df_aggregated_all <- df_all %>%
  group_by(p1_id, p2_id) %>%
  summarise(rg = first(rg)) %>%
  ungroup()

matrix_df_woCorrection <- spread(df_aggregated_all, key = p1_id, value = rg)

row_names <- matrix_df_woCorrection$p2_id
matrix_df_woCorrection$p2_id <- NULL
rownames(matrix_df_woCorrection) <- row_names

matrix_df_woCorrection <- as.matrix(matrix_df_woCorrection)

# Create a matrix containing all adjusted correlation values (with corrected gc and significant p-values)
df_aggregated_all_sign <- df_all_adjusted_significant %>%
  group_by(p1_id, p2_id) %>%
  summarise(rg = first(rg)) %>%
  ungroup()

matrix_df_correction_sign <- spread(df_aggregated_all_sign, key = p1_id, value = rg)

row_names <- matrix_df_correction_sign$p2_id
matrix_df_correction_sign$p2_id <- NULL
rownames(matrix_df_correction_sign) <- row_names
matrix_df_correction_sign[is.na(matrix_df_correction_sign)] <- 0

matrix_df_correction_sign <- as.matrix(matrix_df_correction_sign)

# Save all created df's and matrices 
write.table(df_all, file = paste(output_ld_regression_path, "/complete_df_ldsc.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(df_all_adjusted, file = paste(output_ld_regression_path, "/complete_df_adjusted_ldsc.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(matrix_df, file = paste(output_ld_regression_path, "/matrix_df_ldsc.tsv", sep=""), sep="\t", row.names = TRUE, quote=FALSE)
write.table(matrix_df_woCorrection, file = paste(output_ld_regression_path, "/matrix_df_woCorrection_ldsc.tsv", sep=""), sep="\t", row.names = TRUE, quote=FALSE)
write.table(matrix_df_correction_sign, file = paste(output_ld_regression_path, "/matrix_df_correction_sign_ldsc.tsv", sep=""), sep="\t", row.names=TRUE, quote=FALSE)
write.table(df_all_adjusted_significant, file=paste(output_ld_regression_path, "/df_all_adjusted_significant_ldsc.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
