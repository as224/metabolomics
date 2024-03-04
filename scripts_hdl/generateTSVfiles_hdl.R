# Generates tsv files from the HDL output files to generate heatmaps and networks later on

# Load required libraries
#install.packages("readtext")
library("readtext")
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

output_hdl_path <- "../output_hdl"

pos_inf_counter <- 0
neg_inf_counter <- 0
numeric_counter_inside <- 0
numeric_counter_outside <- 0

# Function that generates tsv from log and adds cols (p1_id, p2_id, id and method)
generate_df <- function(input) {
  log <- readLines(input)
  
  column_names <- c("p1", "p2", "rg", "p", "p1_id", "p2_id")
  
  # Create an empty data frame with known column names
  df <- data.frame(matrix(ncol = length(column_names), nrow = 2))
  colnames(df) <- column_names
  
  for (line in log) {
    # Get gwas1 path
    if (grepl("^gwas1\\.df=", line)) {
      df$p1[1] <- sub("^gwas1\\.df=", "", line)
      df$p1_id[1] <- sub("^.*/output_wrangling/output_", "", gsub(".hdl.rds$", "", df$p1[1]))
    }
    
    # Get gwas2 path
    if (grepl("^gwas2\\.df=", line)) {
      df$p2[1] <- sub("^gwas2\\.df=", "", line)
      df$p2_id[1] <- sub("^.*/output_wrangling/output_", "", gsub(".hdl.rds$", "", df$p2[1]))
    }
    
    # Get genetic covariance 
    if (grepl("^Genetic Correlation:", line)) {
      df$rg[1] <- sub("^Genetic Correlation:\\s+(-?\\d+\\.\\d+).*", "\\1", line)
    }
    
    # Get p value
    if (grepl("^P:", line)) {
      df$p[1] <- sub("^P:", "", line)
    }
    
    #print(df$rg[1])
    
    # Add inverted pair to df (p2/p1)
    df$p1[2] <- df$p2[1]
    df$p2[2] <- df$p1[1]
    df$rg[2] <- df$rg[1]
    df$p[2] <- df$p[1]
    df$p1_id[2] <- df$p2_id[1]
    df$p2_id[2] <- df$p1_id[1]
    
  }
  
  # Add two cols to df
  df$id <- paste(df$p1_id, df$p2_id, sep = "/")
  df$method <- rep("hdl", nrow(df))
  
  # Count calculated rg Infs, -Infs, and numeric values (inside/outside [-1,1])
  if (startsWith(df$rg[1], "Genetic Correlation:  Inf")) {
    pos_inf_counter <<- pos_inf_counter+1 
  }
  else if (startsWith(df$rg[1], "Genetic Correlation:  -Inf")) {
    neg_inf_counter <<- neg_inf_counter+1 
  }
  else if (grepl("^-?\\d+\\.?\\d*$", df$rg[1])) {
    # Inside interval
    if (as.numeric(df$rg[1]) >= -1 && as.numeric(df$rg[1]) <= 1) {
      numeric_counter_inside <<- numeric_counter_inside+1
    }
    # Outside interval
    else {
      numeric_counter_outside <<- numeric_counter_outside+1
    }
  }
  # else {
  #   print(df$rg[1])
  # }
  
  # Replace NA with 0 for genetic correlation 
  df$rg <- as.numeric(df$rg)
  df$rg[is.na(df$rg)] <- 0
  
  # Return the df
  return(df)
}

# Initialize an empty data frame without specifying the number of rows
df_all <- data.frame(matrix(ncol = 10, nrow=0))
colnames(df_all) <- c("p1", "p2", "rg", "se", "p",  "p1_id", "p2_id", "id", "method", "q")

files <- list.files(path = output_hdl_path, pattern = "\\.Rout$", full.names = TRUE)

for (file in files) {
  input_file_path <- file
  output_file_path <- gsub("raw.gwas.Rout", ".tsv", input_file_path) # replace file type
  
  # Call function
  df <- generate_df(input = input_file_path)
  
  df_all <- rbind(df_all, df)
}

# Add diagonal pairs (hdl(gwas1,gwas1)) to df_all
for (element in unique(df_all$p1)) {
  id <- sub("^.*/output_wrangling/output_", "", gsub(".hdl.rds$", "", element))
  id_id <- paste(id, "/", id, sep = "")
  new_row <- c(element, element, 1, NA, id, id, id_id, "hdl", NA)
  
  df_all <- rbind(df_all, new_row)
}

# Adjust the correlation values -> >1=0 and <-1=0 (the highest correlation (diagonal values) is one -> thus a higher value than 1 explains something weird and also <-1)
df_all_adjusted <- df_all
df_all$rg <- as.numeric(df_all$rg) 
df_all_adjusted$rg <- as.numeric(df_all_adjusted$rg)
df_all_adjusted$rg[df_all$rg > 1] <- 0
df_all_adjusted$rg[df_all$rg < -1] <- 0
df_all_adjusted$p <- as.numeric(df_all_adjusted$p)
df_all_adjusted$q <- p.adjust(df_all_adjusted$p, method = "fdr")

# Dataframe with corrected and significant p-values
df_all_adjusted_significant <- df_all_adjusted[df_all_adjusted$p < 0.05 & !(is.na(df_all_adjusted$p)), ]

# Add the diagonal 1 values for left significant gwas
for (element in unique(df_all_adjusted_significant$p1)) {
  id <- sub("^.*/output_wrangling/output_", "", gsub(".hdl.rds$", "", element))
  id_id <- paste(id, "/", id, sep = "")
  new_row <- c(element, element, 1, NA, id, id, id_id, "hdl", NA)
  
  df_all_adjusted_significant <- rbind(df_all_adjusted_significant, new_row)
}

# Create a matrix containing all correlation values (with adjusted correlation values -> >1, <-1 handled)

df_aggregated <- df_all_adjusted %>%
  group_by(p1_id, p2_id) %>%
  summarise(rg = first(rg)) %>%
  ungroup()

matrix_df <- spread(df_aggregated, key = p1_id, value = rg)

row_names <- matrix_df$p2_id
matrix_df$p2_id <- NULL
rownames(matrix_df) <- row_names

matrix_df <- as.matrix(matrix_df)

# Create a matrix containing all correlation values (with real correlation values)

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
matrix_df_correction_sign[] <- lapply(matrix_df_correction_sign, as.numeric) # Convert dataframe to numeric
matrix_df_correction_sign[is.na(matrix_df_correction_sign)] <- 0 # Replace NA values with 0
rownames(matrix_df_correction_sign) <- row_names

matrix_df_correction_sign <- as.matrix(matrix_df_correction_sign)

# Print the amount of Infs, -Infs, genetic correlations outside [-1,1] and inside [-1,1]
print(paste("Amount of GWAS pairs with rg of Inf:", pos_inf_counter))
print(paste("Amount of GWAS pairs with rg of -Inf:", neg_inf_counter))
print(paste("Amount of GWAS pairs with rg outside [-1,1]:", numeric_counter_outside))
print(paste("Amount of GWAS pairs with rg inside [-1,1]:", numeric_counter_inside))

# Save all created df's and matrices 
write.table(df_all, file = paste(output_hdl_path, "/complete_df_hdl.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(df_all_adjusted, file = paste(output_hdl_path, "/complete_df_adjusted_hdl.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(df_all_adjusted_significant, file = paste(output_hdl_path, "/df_all_adjusted_significant_hdl.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
write.table(matrix_df, file = paste(output_hdl_path, "/matrix_df_hdl.tsv", sep=""), sep="\t", row.names = TRUE, quote=FALSE)
write.table(matrix_df_woCorrection, file = paste(output_hdl_path, "/matrix_df_woCorrection_hdl.tsv", sep=""), sep="\t", row.names = TRUE, quote=FALSE)
write.table(matrix_df_correction_sign, file = paste(output_hdl_path, "/matrix_df_correction_sign_hdl.tsv", sep=""), sep="\t", row.names=TRUE, quote=FALSE)
