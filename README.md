# CDCA
A centrality-based community detection method for RNA-Seq data
# Data
Bulk RNA-Seq data for bipolar disorder downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80336, and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53239
Bulk RNA-Seq data for schizophrenia downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138082, and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165604

# Steps
1. Run DE_genes.R to get DE genes for a dataset.
2. Compile homogeneity.R
3. Compile cdca.R
4. Run input.R, which takes the list of DE genes found from Step 1 as input.
