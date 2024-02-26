# This script creates input_rg.txt files (for ldsc and hdl) in order to be able to use the ldsc-network-plot
library(stringr)

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set folder paths
folder_path_ldsc = "../output_ld_regression"
folder_path_hdl = "../output_hdl"

# Read in the significant corrected gc values for ldsc and hdl 
df_ldsc <- read.table(paste(folder_path_ldsc, "/df_all_adjusted_significant_ldsc.tsv", sep=""), header=TRUE)
df_ldsc <- subset(df_ldsc, select = -c(p1 , p2))

df_hdl <- read.table(paste(folder_path_hdl, "/df_all_adjusted_significant_hdl.tsv", sep=""), header=TRUE)
df_hdl <- subset(df_hdl, select = -c(p1 , p2))

# Adjust column names 
colnames(df_ldsc)[colnames(df_ldsc) == "p1_id"] <- "p1"
colnames(df_ldsc)[colnames(df_ldsc) == "p2_id"] <- "p2"

colnames(df_hdl)[colnames(df_hdl) == "p1_id"] <- "p1"
colnames(df_hdl)[colnames(df_hdl) == "p2_id"] <- "p2"

trait_category_mapping <- c("IL" = "Cytokine",
                            "MS" = "Autoimmune",
                            "ALS" = "Autoimmune",
                            "bNGF" = "Protein",
                            "CTACK" = "Cytokine",
                            "eotaxin" = "Cytokine",
                            "GCSF" = "Protein",
                            "FGF" = "Protein",
                            "GROa" = "Cytokine",
                            "IFNg" = "Cytokine",
                            "ip10" = "Protein",
                            "MCP" = "Protein",
                            "MCSF" = "Cytokine",
                            "MIF" = "Cytokine",
                            "MIG" = "Cytokine",
                            "MIP" = "Cytokine",
                            "PDGFbb" = "Protein",
                            "RANTES" = "Cytokine",
                            "SCF" = "Cytokine",
                            "SCGFb" = "Protein",
                            "SDF1" = "Cytokine",
                            "TNF" = "Cytokine",
                            "TRAIL" = "Cytokine",
                            "T2D" = "Autoimmune",
                            "PD" = "Autoimmune",
                            "RA" = "Inflammatory",
                            "DLB" = "Inflammatory",
                            "SCZ" = "Inflammatory",
                            "C3" = "Protein",
                            "HGF" = "Protein",
                            "MCP1" = "Protein",
                            "CRP" = "Protein",
                            "APH"  = "Autoimmune",
                            "tryptophan" = "Aminoacid",
                            "histidine" = "Aminoacid",
                            "phenylalanine" = "Aminoacid",
                            "isoleucine" = "Aminoacid",
                            "threonine" = "Aminoacid",
                            "tyrosine" = "Aminoacid",
                            "arginine" = "Aminoacid",
                            "proline" = "Aminoacid",
                            "acetylcarnitine" = "Other",
                            "glutamate" = "Aminoacid",
                            "alanine" = "Aminoacid",
                            "asparagine" = "Aminoacid",
                            "cholesterol" = "Other",
                            "pyruvate" = "Other",
                            "citrate" = "Other",
                            "urea" = "Other",
                            "glycerol" = "Other")

# Prerequisite for mapping the categories to the traits 
df_ldsc$p1_helper <- df_ldsc$p1
df_ldsc$p2_helper <- df_ldsc$p2

df_ldsc$p1_helper <- sapply(str_split(df_ldsc$p1_helper, "_"), function(x) x[[2]])
df_ldsc$p2_helper <- sapply(str_split(df_ldsc$p2_helper, "_"), function(x) x[[2]])

df_hdl$p1_helper <- df_hdl$p1
df_hdl$p2_helper <- df_hdl$p2

df_hdl$p1_helper <- sapply(str_split(df_hdl$p1_helper, "_"), function(x) x[[2]])
df_hdl$p2_helper <- sapply(str_split(df_hdl$p2_helper, "_"), function(x) x[[2]])


for (key in names(trait_category_mapping)) {
  
  # Map ldsc categories and traits
  indices_p1 <- grepl(paste0("^", key), df_ldsc$p1_helper)
  indices_p2 <- grepl(paste0("^", key), df_ldsc$p2_helper)
  
  df_ldsc$p1_category[indices_p1] <- trait_category_mapping[key]
  df_ldsc$p2_category[indices_p2] <- trait_category_mapping[key]
  
  # Map hdl categories and traits
  indices_p1 <- grepl(paste0("^", key), df_hdl$p1_helper)
  indices_p2 <- grepl(paste0("^", key), df_hdl$p2_helper)
  
  df_hdl$p1_category[indices_p1] <- trait_category_mapping[key]
  df_hdl$p2_category[indices_p2] <- trait_category_mapping[key]
}

df_ldsc <- subset(df_ldsc, p1 != p2)
df_ldsc <- subset(df_ldsc, select = c("p1_category", "p1", "p2_category", "p2", "rg", "se", "z", "p", "q"))

df_hdl <- subset(df_hdl, p1 != p2)
df_hdl <- subset(df_hdl, select = c("p1_category", "p1", "p2_category", "p2", "rg", "se", "z", "p", "q"))

# Exclude GWAS pairs where gc=0 
# -> Either it was >1 or <-1 or it was NA or truly totally independent 
# --> no handling planned 
# also if truly independent (0) then it would just be an independent point in the plot))


df_ldsc <- df_ldsc[df_ldsc$rg != 0, ] 
df_hdl <- df_hdl[df_hdl$rg != 0, ] 

# Save as txt files 
write.table(df_ldsc, file = paste(folder_path_ldsc, "/input_rg_ldsc.txt", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df_hdl, file = paste(folder_path_hdl, "/input_rg_hdl.txt", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


