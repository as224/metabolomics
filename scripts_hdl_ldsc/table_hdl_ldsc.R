# The goal of this script is to create two tables. Once including only significant gc values according to ldsc and once according to hdl. 
library(grid)
library(gridExtra)

# Set the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set paths to df_all_adjusted_significant_ldsc.tsv and df_all_adjusted_significant_hdl.tsv 
# These files contain the corrected gc values (gc>1 -> 0, gc<-1 -> 0). Also they contain a fdr value and only significant p values (p<0.05).
input_file_ldsc_sign <- "../output_ld_regression/df_all_adjusted_significant_ldsc.tsv"
input_file_hdl_sign <- "../output_hdl/df_all_adjusted_significant_hdl.tsv"

# Set paths to complete_df_adjusted_ldsc.tsv and complete_df_adjusted_hdl.tsv
# These files contain all corrected gc values (no filter for significance).
input_file_ldsc_adjusted <- "../output_ld_regression/complete_df_adjusted_ldsc.tsv"
input_file_hdl_adjusted <- "../output_hdl/complete_df_adjusted_hdl.tsv"

# Read all files in and transform the data frame the way we need it. 
df_ldsc_sign <- read.table(input_file_ldsc_sign, header=TRUE)
df_ldsc_sign <- subset(df_ldsc_sign, p1 != p2)
df_ldsc_sign <- subset(df_ldsc_sign, select = c("p1_id", "p2_id", "rg", "p"))
colnames(df_ldsc_sign) <- c("GWAS 1", "GWAS 2", "gc (LDSC)", "p-value (LDSC)")

df_ldsc_adj <- read.table(input_file_ldsc_adjusted, header=TRUE)
df_ldsc_adj <- subset(df_ldsc_adj, p1 != p2)
df_ldsc_adj <- subset(df_ldsc_adj, select = c("p1_id", "p2_id", "rg", "p"))
colnames(df_ldsc_adj) <- c("GWAS 1", "GWAS 2", "gc (LDSC)", "p-value (LDSC)")


df_hdl_sign <- read.table(input_file_hdl_sign, header=TRUE)
df_hdl_sign <- subset(df_hdl_sign, p1 != p2)
df_hdl_sign <- subset(df_hdl_sign, select = c("p1_id", "p2_id", "rg", "p"))
colnames(df_hdl_sign) <- c("GWAS 1", "GWAS 2", "gc (HDL)", "p-value (HDL)")

df_hdl_adj <- read.table(input_file_hdl_adjusted, header=TRUE)
df_hdl_adj <- subset(df_hdl_adj, p1 != p2)
df_hdl_adj <- subset(df_hdl_adj, select = c("p1_id", "p2_id", "rg", "p"))
colnames(df_hdl_adj) <- c("GWAS 1", "GWAS 2", "gc (HDL)", "p-value (HDL)")


# Merge the data frames
# Only significant gc values according to LDSC + comparisson to gc values of HDL (no significance necessary)
ldsc_sign_hdl_adj <- merge(df_ldsc_sign, df_hdl_adj, by = c("GWAS 1", "GWAS 2"), all.x = TRUE)
# Only significant gc values according to HDL + comparisson to gc values of LDSC (no significance necessary) 
hdl_sign_ldsc_adj <- merge(df_hdl_sign, df_ldsc_adj, by = c("GWAS 1", "GWAS 2"), all.x = TRUE)

ldsc_sign_hdl_adj[is.na(ldsc_sign_hdl_adj)] <- 0  # Replace NA with 0 or any other desired value
hdl_sign_ldsc_adj[is.na(hdl_sign_ldsc_adj)] <- 0  # Replace NA with 0 or any other desired value

ldsc_sign_hdl_adj <- subset(ldsc_sign_hdl_adj, select = c("GWAS 1", "GWAS 2", "gc (HDL)", "gc (LDSC)", "p-value (HDL)", "p-value (LDSC)"))
hdl_sign_ldsc_adj <- subset(hdl_sign_ldsc_adj, select = c("GWAS 1", "GWAS 2", "gc (HDL)", "gc (LDSC)", "p-value (HDL)", "p-value (LDSC)"))

ldsc_sign_hdl_adj <- ldsc_sign_hdl_adj[ldsc_sign_hdl_adj$`p-value (LDSC)` != 0, ]
hdl_sign_ldsc_adj <- hdl_sign_ldsc_adj[hdl_sign_ldsc_adj$`p-value (HDL)` != 0, ]

resfactor=20
png("../plots/table_ldsc_sign_hdl_adj_head.png", res = 72*resfactor, height=500*resfactor, width=850*resfactor) 
grid.table(head(ldsc_sign_hdl_adj, 20))
dev.off()

png("../plots/table_hdl_sign_ldsc_adj_head.png", res = 72*resfactor, height=80*resfactor, width=850*resfactor)
grid.table(head(hdl_sign_ldsc_adj, 10))
dev.off()



