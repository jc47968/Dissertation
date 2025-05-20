library(dplyr)                                                 
library(readr)
library(data.table)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\03 Merged Data]")
n_data <- fread("phenotype_data.csv")
corr_data <- fread("correlation_results_wide.csv")

# Ensure both datasets are data.tables
setDT(n_data)
setDT(corr_data)


dim(n_data) #[1] 568  15
dim(corr_data) #[1]  1053 76638


# Merge datasets on SUB_ID (n_data) and Subject (corr_data)
merged_data <- merge(
  n_data,
  corr_data,
  by.x = "SUB_ID",
  by.y = "Subject",
  all = FALSE #inner join
)

# Remove the Subject column if still present (it's merged as SUB_ID)
merged_data[, Subject := NULL]

dim(merged_data)#[1]   486 76652

write.csv(merged_data, "merged_data.csv")
