# Run this file after compiling homogeneity.R and cdca.R

setwd("DE_genes/")

#correlation matrix
# SZ
data_D1_gse165604 <- read.csv(file = "SZ/corr_D1_gse165604.csv", stringsAsFactors = FALSE)
dim(data_D1_gse165604)  # 146 x 147

community_component <- eigen_cen_community(data_D1_gse165604, "corr_D1_gse165604")



data_C1_gse165604 <- read.csv(file = "SZ/corr_C1_gse165604.csv", stringsAsFactors = FALSE)
dim(data_C1_gse165604)  # 146 x 147

community_component2 <- eigen_cen_community(data_C1_gse165604, "corr_C1_gse165604")

#----

data_D1_gse138082 <- read.csv(file = "SZ/corr_D1_gse138082.csv", stringsAsFactors = FALSE)
dim(data_D1_gse138082)  # 225 226

community_component3 <- eigen_cen_community(data_D1_gse138082, "corr_D1_gse138082")



data_C1_gse138082 <- read.csv(file = "SZ/corr_C1_gse138082.csv", stringsAsFactors = FALSE)
dim(data_C1_gse138082)  #  225 226

community_component4 <- eigen_cen_community(data_C1_gse138082, "corr_C1_gse138082")

#----

# BD
data_D1_gse53239 <- read.csv(file = "BD/corr_D1_gse53239.csv", stringsAsFactors = FALSE)
dim(data_D1_gse53239)  # 37 38

community_component5 <- eigen_cen_community(data_D1_gse53239, "corr_D1_gse53239")



data_C1_gse53239 <- read.csv(file = "BD/corr_C1_gse53239.csv", stringsAsFactors = FALSE)
dim(data_C1_gse53239)  # 37 38

community_component6 <- eigen_cen_community(data_C1_gse53239, "corr_C1_gse53239")

#----


data_D1_gse80336 <- read.csv(file = "BD/corr_D1_gse80336.csv", stringsAsFactors = FALSE)
dim(data_D1_gse80336)  # 47 48

community_component7 <- eigen_cen_community(data_D1_gse80336, "corr_D1_gse80336")



data_C1_gse80336 <- read.csv(file = "BD/corr_C1_gse80336.csv", stringsAsFactors = FALSE)
dim(data_C1_gse80336)  # 47 48

community_component8 <- eigen_cen_community(data_C1_gse80336, "corr_C1_gse80336")


