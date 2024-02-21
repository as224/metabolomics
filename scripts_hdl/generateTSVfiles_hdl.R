# Load required libraries
#install.packages("readtext")
library("readtext")
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#output_hdl_path <- "../output_hdl"
output_hdl_path <- "../large_gwas/output_hdl"

# Function that generates tsv from log and adds cols (p1_id, p2_id, id and method)
#' Title
#'
#' @param input 
#'
#' @return
#' @export
#'
#' @examples
generate_df <- function(input) {
  log <- readLines(input)
  
  column_names <- c("p1", "p2", "rg", "p", "p1_id", "p2_id")
  
  # Create an empty data frame with known column names
  df <- data.frame(matrix(ncol = length(column_names), nrow = 1))
  colnames(df) <- column_names
  
  for (line in log) {
    # Get gwas1 path
    if (grepl("^gwas1\\.df=", line)) {
      # cols$p1 <- sub("^gwas1\\.df=", "../", line)
      # cols$p1_id <- sub("^.*?/pmid[0-9]+_(.*?)\\.hdl\\.rds$", "\\1", cols$p1)
      df$p1[1] <- sub("^gwas1\\.df=", "../", line)
      df$p1_id[1] <- sub("^.*?/pmid[0-9]+_(.*?)\\.hdl\\.rds$", "\\1", df$p1[1])
    }
    
    # Get gwas2 path
    if (grepl("^gwas2\\.df=", line)) {
      #cols$p2 <- sub("^gwas2\\.df=", "../", line)
      #cols$p2_id <- sub("^.*?/pmid[0-9]+_(.*?)\\.hdl\\.rds$", "\\1", cols$p2)
      df$p2[1] <- sub("^gwas2\\.df=", "../", line)
      df$p2_id[1] <- sub("^.*?/pmid[0-9]+_(.*?)\\.hdl\\.rds$", "\\1", df$p2[1])
    }
    
    # Get genetic covariance 
    if (grepl("^Genetic Correlation:", line)) {
      #cols$rg <- sub("^Genetic Covariance:", "", line)
      df$rg <- sub("^Genetic Correlation:\\s+(-?\\d+\\.\\d+).*", "\\1", line)
    }
    
    # Überprüfen, ob die Zeile Informationen zum p-Wert enthält
    if (grepl("^P:", line)) {
      #cols$p <- sub("^P:", "", line)
      df$p <- sub("^P:", "", line)
    }
  }

  # Umwandeln der Liste in einen Datenrahmen
  #df <- as.data.frame(t(cols))
  
  # Add two cols to df
  df$id <- paste(df$p1_id, df$p2_id, sep = "/")
  df$method <- rep("hdl", nrow(df))
  
  # Replace NA with 0 for genetic correlation 
  df$rg <- as.numeric(df$rg)
  df$rg[is.na(df$rg)] <- 0
  
  # Add FDR 
  df$q <- p.adjust(df$p, method = "fdr")
  
  # Return the df
  return(df)
}

# Initialize an empty data frame without specifying the number of rows
df_all <- data.frame(matrix(ncol = 9, nrow=0))
colnames(df_all) <- c("p1", "p2", "rg", "se", "p",  "p1_id", "p2_id", "id", "method")

files <- list.files(path = output_hdl_path, pattern = "\\.Rout$", full.names = TRUE)

for (file in files) {
  input_file_path <- file
  output_file_path <- gsub("raw.gwas.Rout", ".tsv", input_file_path) # replace file type
  
  # Call function
  df <- generate_df(input = input_file_path)
  
  df_all <- rbind(df_all, df)
}

print(df_all)
class(df_all$p1)

# Adjust the correlation values -> >1=0 and <-1=0 (the highest correlation (diagnoal values) is one -> thus a higher value than 1 explains something weird and also <-1)
df_all_adjusted <- df_all  
df_all_adjusted$rg[df_all$rg > 1] <- 0
df_all_adjusted$rg[df_all$rg < -1] <- 0

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

# Save all created df's and matrices 
write.table(df_all, file = paste(output_hdl_path, "/complete_df.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(df_all_adjusted, file = paste(output_hdl_path, "/complete_df_adjusted.tsv", sep="") , sep = "\t", row.names = FALSE, quote=FALSE)
write.table(matrix_df, file = paste(output_hdl_path, "/matrix_df.tsv", sep=""), row.names = TRUE, quote=FALSE)
write.table(matrix_df_woCorrection, file = paste(output_hdl_path, "/matrix_df_woCorrection.tsv", sep=""), row.names = TRUE, quote=FALSE)
