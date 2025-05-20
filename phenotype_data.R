library(dplyr)                                                 
library(readr)
library(data.table)
library(reshape2)
library(ggplot2)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\01 Phenotype Data]")
p_data <- fread("Phenotypic_V1_0b_v1.csv")
print(dim(p_data)) #[1] 1112   74

# Select the needed columns
p_data1 <- p_data %>% dplyr::select(SITE_ID, SUB_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN,
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
  x <- as.numeric(as.character(x))
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(x)
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

##
# Create a comparison table with summary statistics before and after normalization
comparison_table <- data.frame(Variable = columns_to_check)

comparison_table <- comparison_table %>%
  mutate(
    Mean_Before = sapply(columns_to_check, function(col) mean(p_data1[[col]], na.rm = TRUE)),
    Min_Before = sapply(columns_to_check, function(col) min(p_data1[[col]], na.rm = TRUE)),
    Max_Before = sapply(columns_to_check, function(col) max(p_data1[[col]], na.rm = TRUE)),
    Mean_After = sapply(columns_to_check, function(col) mean(n_data[[col]], na.rm = TRUE)),
    Min_After = sapply(columns_to_check, function(col) min(n_data[[col]], na.rm = TRUE)),
    Max_After = sapply(columns_to_check, function(col) max(n_data[[col]], na.rm = TRUE))
  )

# Print the table
print(comparison_table)

# Save the table as a CSV file
write.csv(comparison_table, "Normalization_Comparison.csv", row.names = FALSE)

##
# Melt the data for visualization
p_data1_melted <- melt(p_data1[, ..columns_to_check], variable.name = "Variable", value.name = "Original_Value")
p_data1_melted$Normalization <- "Before"

n_data_melted <- melt(n_data[, ..columns_to_check], variable.name = "Variable", value.name = "Normalized_Value")
n_data_melted$Normalization <- "After"

# Merge before and after data
comparison_data <- merge(p_data1_melted, n_data_melted, by = c("Variable", "Normalization"))

# Create the dual-axis plot
ggplot() +
  # Boxplot for original values (Left Y-Axis)
  geom_boxplot(data = p_data1_melted, aes(x = Variable, y = Original_Value, fill = Normalization), alpha = 0.5) +
  
  # Boxplot for normalized values (Right Y-Axis)
  geom_boxplot(data = n_data_melted, aes(x = Variable, y = Normalized_Value * max(p_data1_melted$Original_Value, na.rm = TRUE), fill = Normalization), alpha = 0.5) +
  
  scale_y_continuous(
    name = "Original Scale",
    sec.axis = sec_axis(~ . / max(p_data1_melted$Original_Value, na.rm = TRUE), name = "Normalized Scale (0-1.00)")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Comparison of Variables Before and After Normalization",
       x = "Variable",
       fill = "Normalization",
       color = "Normalization") +
  scale_fill_manual(values = c("Before" = "dark gray", "After" = "light gray"))

n_data$DX_GROUP <- as.factor(n_data$DX_GROUP)

# Boxplot: AVG_SOCIAL and AVG_COMM by DX_GROUP
ggplot(n_data, aes(x = DX_GROUP)) +
  geom_boxplot(aes(y = AVG_SOCIAL, fill = "AVG_SOCIAL"), alpha = 0.6) +
  geom_boxplot(aes(y = AVG_COMM, fill = "AVG_COMM"), alpha = 0.6) +
  labs(title = "Comparison of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "DX_GROUP",
       y = "Score",
       fill = "Variable") +
  theme_minimal() +
  scale_fill_manual(values = c("AVG_SOCIAL" = "dark gray", "AVG_COMM" = "light gray")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reshape data for ggplot
data_long <- n_data %>%
  pivot_longer(cols = c(AVG_SOCIAL, AVG_COMM), names_to = "Variable", values_to = "Score")

# Create density plot with proper legend for line types
ggplot(data_long, aes(x = Score, color = DX_GROUP, linetype = Variable)) +
  geom_density(size = 1.2) +
  scale_linetype_manual(values = c("AVG_SOCIAL" = "solid", "AVG_COMM" = "dotted")) +
  labs(title = "Density Plot of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "Score",
       y = "Density",
       color = "DX_GROUP",
       linetype = "Variable") +  
  theme_minimal() +
  theme(legend.position = "top")# Density plot: AVG_SOCIAL and AVG_COMM by DX_GROUP
ggplot() +
  geom_density(data = n_data, aes(x = AVG_SOCIAL, color = DX_GROUP, linetype = "AVG_SOCIAL"), size = 1.2) +
  geom_density(data = n_data, aes(x = AVG_COMM, color = DX_GROUP, linetype = "AVG_COMM"), size = 1.2, linetype = "dotted") +
  scale_linetype_manual(values = c("AVG_SOCIAL" = "solid", "AVG_COMM" = "dotted"),
                        name = "Variable", labels = c("AVG_SOCIAL (Solid)", "AVG_COMM (Dotted)")) +
  labs(title = "Density Plot of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "Score",
       y = "Density",
       color = "DX_GROUP") +
  theme_minimal() +
  theme(legend.position = "top")

