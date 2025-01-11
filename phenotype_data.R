library(dplyr)                                                 
library(readr)
library(data.table)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\01 Phenotype Data]")
p_data <- fread("Phenotypic_V1_0b_v1.csv")
print(dim(p_data)) #[1] 1112   74

# Select the needed columns
p_data1 <- p_data %>% select(SITE_ID, SUB_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN,
                            SEX, ADI_R_SOCIAL_TOTAL_A, ADI_R_VERBAL_TOTAL_BV,
                            ADOS_COMM, ADOS_SOCIAL, SRS_COGNITION, SRS_COMMUNICATION)
print(dim(p_data1))#[1] 1112   12

# Get the column range between "ADI_R_SOCIAL_TOTAL_A" and "SRS_COMMUNICATION"
start_col <- which(names(p_data1) == "ADI_R_SOCIAL_TOTAL_A")
end_col <- which(names(p_data1) == "SRS_COMMUNICATION")

# Subset columns in the range
columns_to_check <- names(p_data1)[start_col:end_col]

# Replace all -9999 values with NA in the entire data frame
p_data1[p_data1 == -9999] <- NA

# Remove rows where all data in the specified columns are blank or NA
cleaned_data <- p_data1 %>%
  filter(
    rowSums(
      sapply(select(., all_of(columns_to_check)), function(x) !is.na(x) & x != "")
    ) > 0
  )
print(dim(cleaned_data)) #[1] 568  12

# Normalize the data in the selected columns
setDT(cleaned_data)
n_data <- cleaned_data
n_data[, (columns_to_check) := lapply(.SD, function(x) {
  # Convert to numeric if necessary
  x <- as.numeric(as.character(x))
  # Avoid issues with NA or constant values
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(x) # Return original column if all values are NA or constant
  }
  # Normalize (scale to 0-1)
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}), .SDcols = columns_to_check]

# Create the AVG_SOCIAL column
n_data[, AVG_SOCIAL := rowMeans(.SD, na.rm = TRUE), 
     .SDcols = c("ADI_R_SOCIAL_TOTAL_A", "ADOS_SOCIAL", "SRS_COGNITION")]

# Create the AVG_COMM column
n_data[, AVG_COMM := rowMeans(.SD, na.rm = TRUE), 
     .SDcols = c("ADI_R_VERBAL_TOTAL_BV", "ADOS_COMM", "SRS_COMMUNICATION")]

write.csv(n_data, "phenotype_data.csv")
