# This script creates the traitlist.txt file as input for ldsc-network-plot (for ldsc and hdl)
library(stringr)

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# The traitlist for hdl and ldsc is the same. Here I take the unique IDs from ldsc output but one can also take that from hdl output.
folder_path = "../output_ld_regression"
df_pre_trait <- read.table(paste(folder_path, "/complete_df_adjusted_ldsc.tsv", sep=""), header=TRUE)

# Get unique values from df_pre_trait$p1_id
unique_p1_id <- unique(df_pre_trait$p1_id)

# Create an empty data frame
df_trait <- data.frame(
  CATEGORY = numeric(length(unique_p1_id)),
  TRAIT = character(length(unique_p1_id)),
  COLOR = logical(length(unique_p1_id))
)

# Assign trait values
df_trait$TRAIT <- unique_p1_id

df_trait$trait_helper <- df_trait$TRAIT
df_trait$trait_helper <- sapply(str_split(df_trait$trait_helper, "_"), function(x) x[[2]])

# Create a lookup table for mapping TRAIT to CATEGORY
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
                            "glycerol" = "Other",
                            "T2d" = "Autoimmune",
                            "LBD" = "Inflammatory")

for (key in names(trait_category_mapping)) {
  # Find indices where TRAIT starts with the key
  indices <- grepl(paste0("^", key), df_trait$trait_helper)
  
  # Update CATEGORY with corresponding category from trait_category_mapping
  df_trait$CATEGORY[indices] <- trait_category_mapping[key]
}



df_trait$COLOR <- ifelse(startsWith(df_trait$CATEGORY, "Cytokine"), "red", 
                         ifelse(startsWith(df_trait$CATEGORY, "Protein"), "lightgreen", 
                                ifelse(startsWith(df_trait$CATEGORY, "Autoimmune"), "lightblue", 
                                       ifelse(startsWith(df_trait$CATEGORY, "Inflammatory"), "violet",
                                              ifelse(startsWith(df_trait$CATEGORY, "Aminoacid"), "lightpink",
                                                     "white")))))

df_trait <- subset(df_trait, select = c("CATEGORY", "TRAIT", "COLOR"))


write.table(df_trait, file = paste(folder_path, "/traitlist.txt", sep="") ,sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



