# This script creates a input_rg.txt file in order to be able to use ldsc-network-plot
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder_path = "../output_ld_regression"

df_adjusted <- read.table(paste(folder_path, "/complete_df_adjusted.tsv", sep=""), header=TRUE)
df_adjusted <- subset(df_adjusted, select = -c(p1 , p2))

# Adjust column names 
colnames(df_adjusted)[colnames(df_adjusted) == "p1_id"] <- "p1"
colnames(df_adjusted)[colnames(df_adjusted) == "p2_id"] <- "p2"

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
                            "CRP" = "Protein")

# map the categories to the traits 
for (key in names(trait_category_mapping)) {
  indices_p1 <- grepl(paste0("^", key), df_adjusted$p1)
  indices_p2 <- grepl(paste0("^", key), df_adjusted$p2)
  
  df_adjusted$p1_category[indices_p1] <- trait_category_mapping[key]
  df_adjusted$p2_category[indices_p2] <- trait_category_mapping[key]
}

df_adjusted <- subset(df_adjusted, p1 != p2)
df_adjusted <- subset(df_adjusted, select = c("p1_category", "p1", "p2_category", "p2", "rg", "se", "z", "p", "q"))

df_filtered <- df_adjusted[df_adjusted$rg != 0, ] # exclude gwas-gwas pair where gc=0 (either it was >1 or <-1 or it was NA or truly totally independent --> no handling planned (also if truly independen (0) then it would just be an independent point in the plot))

write.table(df_filtered, file = paste(folder_path, "/input_rg_categories_filtered.txt", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


