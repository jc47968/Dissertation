#Permutation MANCOVA for Communication
library(dplyr)                                                 
library(data.table)
library(car)
library(vegan)  # For permutation tests

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\06.2 Perm MANCOVA_Alldata]")

# Load the dataset
vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

# Define the dependent variables (PC1 to PC20)
# Define dependent variables
dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, ..dv_names]
dependent_vars_v <- vdata[, ..dv_names]

set.seed(123) 
perm_mancova_p <- anova.cca(capscale(as.matrix(dependent_vars_p) ~ DX_GROUP * AVG_COMM_INT +
                                       AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                                       SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                                     data = pdata), by = "margin", permutations = 10000)
print(perm_mancova_p)
write.csv(perm_mancova_p, "PermMANCOVA_alldata_promax_comm.csv")

perm_mancova_v <- anova.cca(capscale(as.matrix(dependent_vars_v) ~ DX_GROUP * AVG_COMM_INT +
                                       AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                                       SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                                     data = vdata), by = "margin", permutations = 10000)
print(perm_mancova_v)
write.csv(perm_mancova_v, "PermMANCOVA_alldata_varimax_comm.csv")


