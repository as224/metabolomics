# Metabolomics

The main goal of this profect is exploring genetic correlations and molecular overlaps in inflammatory and autoimmune diseases to prioritize molecular markers and pathways. 


## LDSC 
The LD Score Regression method, designed to estimate heritability and genetic correlation between different GWAS summary statistics, is implemented as a tool. The main repository for it can be found here: https://github.com/bulik/ldsc.
In our analysis, we applied their algorithm to 67 GWAS studies. The list of their ID's can be found at /input_files/GWAS_ids.tsv. To maintain the proper folder structure in your cloned GitHub repository, please clone their GitHub into our main repository with the following structure: metabolomics/ldsc.

Below one can familiarize with our workflow for utilizing LDSC and HDL. It is crucial to follow these steps in order to reproduce our results.


<iframe src="https://docs.google.com/viewer?url=https://github.com/as224/metabolomics/blob/main/workflow_ldsc_hdl.pdf&embedded=true" width="100%" height="600px"></iframe>
