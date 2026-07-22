## ================================================================================
## CONFIG — paths you may need to update
## ================================================================================

# This folder: all pipeline inputs/outputs live here unless noted otherwise
setwd(r"[C:\Users\jason\OneDrive\NCU\Dissertation Dataset\Disseration Analysis\Full_Stack - Copy]")

# Raw per-subject CC400 connectivity time-series CSVs (only needed for Stage 2)
RAW_CC400_DIR <- r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\CC400_CPAC]"

# Combined CC400 time series file that Stage 2 produces / Stage 3 consumes.
# NOTE: this already exists on disk from a prior run (~840 MB), so Stage 2
# below skips re-combining 1053 files unless you delete it.
CC400_COMBINED_FILE <- r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\CC400_combined.csv]"

# CC400 ROI label lookup (only needed for Stages 9-12 brain-network visuals)
ROI_LABELS_FILE <- r"[C:\Users\jason\OneDrive\NCU\Dissertation Dataset\BrainAtlasLabels\CC400_ROI_labels.csv]"

## ================================================================================
## LIBRARIES — consolidated from all 27 scripts (deduplicated)
## ================================================================================

library(dplyr)
library(data.table)
library(readr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(corrplot)
library(GPArotation)
library(psych)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggimage)
library(tidyverse)
library(MVN)
library(car)
library(biotools)
library(vegan)
library(e1071)
library(caret)
library(pROC)
library(multcomp)
library(rstatix)

run_stage <- function(n, title, expr) {
  message(sprintf("\n===== STAGE %d: %s =====", n, title))
  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) message(sprintf("STAGE %d FAILED: %s", n, conditionMessage(e)))
  )
}

## ================================================================================
## STAGE 1 — Phenotype Data Preparation  (source: phenotype_data.R)
## Input: Phenotypic_V1_0b_v1.csv  ->  Output: phenotype_data.csv
## ================================================================================
run_stage(1, "Phenotype Data Preparation", {

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

print(comparison_table)
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
print(ggplot() +
  geom_boxplot(data = p_data1_melted, aes(x = Variable, y = Original_Value, fill = Normalization), alpha = 0.5) +
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
  scale_fill_manual(values = c("Before" = "dark gray", "After" = "light gray")))

n_data$DX_GROUP <- as.factor(n_data$DX_GROUP)

# Boxplot: AVG_SOCIAL and AVG_COMM by DX_GROUP
print(ggplot(n_data, aes(x = DX_GROUP)) +
  geom_boxplot(aes(y = AVG_SOCIAL, fill = "AVG_SOCIAL"), alpha = 0.6) +
  geom_boxplot(aes(y = AVG_COMM, fill = "AVG_COMM"), alpha = 0.6) +
  labs(title = "Comparison of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "DX_GROUP",
       y = "Score",
       fill = "Variable") +
  theme_minimal() +
  scale_fill_manual(values = c("AVG_SOCIAL" = "dark gray", "AVG_COMM" = "light gray")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Reshape data for ggplot
data_long <- n_data %>%
  pivot_longer(cols = c(AVG_SOCIAL, AVG_COMM), names_to = "Variable", values_to = "Score")

# Density plot: AVG_SOCIAL and AVG_COMM by DX_GROUP
print(ggplot(data_long, aes(x = Score, color = DX_GROUP, linetype = Variable)) +
  geom_density(size = 1.2) +
  scale_linetype_manual(values = c("AVG_SOCIAL" = "solid", "AVG_COMM" = "dotted")) +
  labs(title = "Density Plot of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "Score",
       y = "Density",
       color = "DX_GROUP",
       linetype = "Variable") +
  theme_minimal() +
  theme(legend.position = "top"))

print(ggplot() +
  geom_density(data = n_data, aes(x = AVG_SOCIAL, color = DX_GROUP, linetype = "AVG_SOCIAL"), size = 1.2) +
  geom_density(data = n_data, aes(x = AVG_COMM, color = DX_GROUP, linetype = "AVG_COMM"), size = 1.2, linetype = "dotted") +
  scale_linetype_manual(values = c("AVG_SOCIAL" = "solid", "AVG_COMM" = "dotted"),
                        name = "Variable", labels = c("AVG_SOCIAL (Solid)", "AVG_COMM (Dotted)")) +
  labs(title = "Density Plot of AVG_SOCIAL and AVG_COMM by DX_GROUP",
       x = "Score",
       y = "Density",
       color = "DX_GROUP") +
  theme_minimal() +
  theme(legend.position = "top"))

})

## ================================================================================
## STAGE 2 — CC400 Time-Series Combination  (source: CPAC_400_Combine.R)
## Input: raw per-subject CSVs in RAW_CC400_DIR  ->  Output: CC400_COMBINED_FILE
## Skipped automatically if CC400_COMBINED_FILE already exists (it does, ~840MB).
## ================================================================================
run_stage(2, "CC400 Time-Series Combination", {

if (file.exists(CC400_COMBINED_FILE)) {
  message("CC400_combined.csv already exists at ", CC400_COMBINED_FILE, " — skipping combine step.")
} else {

  path1 <- RAW_CC400_DIR
  files <- list.files(pattern = ".csv", path = path1, full.names = TRUE)
  length(files) #1053

  # Function to read multiple CSV files and replace the first column with the filename
  read_and_replace_first_column <- function(file_path) {
    filename <- basename(file_path)
    data <- read_csv(file_path, col_types = cols(.default = "c"))  # Read all columns as character type
    data[[1]] <- filename
    return(data)
  }

  csv_files <- list.files(path1, pattern = "*.csv", full.names = TRUE)
  all_data <- bind_rows(lapply(csv_files, read_and_replace_first_column))
  print(all_data)

  write_csv(all_data, CC400_COMBINED_FILE)
}

})

## ================================================================================
## STAGE 3 — Functional Connectivity Correlation Analysis  (source: Correlation_Analysis.R)
## Input: CC400_COMBINED_FILE  ->  Output: correlation_results_wide.csv
## ================================================================================
run_stage(3, "Functional Connectivity Correlation Analysis", {

data <- read.csv(CC400_COMBINED_FILE)

# Example correlation function to apply for each group
correlation_analysis <- function(df) {
  data_for_correlation <- df[, -1] %>%
    select_if(is.numeric)

  if (ncol(data_for_correlation) > 1) {
    cor_matrix <- cor(data_for_correlation, use = "pairwise.complete.obs")
    cor_matrix[upper.tri(cor_matrix)] <- NA
    cor_matrix <- as.data.frame(as.table(cor_matrix))
    cor_matrix <- cor_matrix[!is.na(cor_matrix$Freq), ]
    cor_matrix <- cor_matrix[cor_matrix$Var1 != cor_matrix$Var2, ]
    colnames(cor_matrix) <- c("Variable1", "Variable2", "Correlation")
  } else {
    cor_matrix <- NA
  }
  return(cor_matrix)
}

# Group by 'Site' and 'Subject' fields and apply the correlation function to each group
correlation_results <- data %>%
  group_by(Site, Subject) %>%
  do(correlation_analysis = correlation_analysis(.)) %>%
  unnest(cols = c(correlation_analysis))

correlation_results <- correlation_results %>%
  mutate(Variable_Pair = paste(Variable1, Variable2, sep = "_")) %>%
  select(-Variable1, -Variable2)

correlation_results_wide <- correlation_results %>%
  spread(key = Variable_Pair, value = Correlation)

corr_data_noNA <- na.omit(correlation_results_wide) #remove rows with NA

write.csv(corr_data_noNA, "correlation_results_wide.csv", row.names = FALSE)

##Illustration - Correlation Matrix of First Observation
first_observation <- data %>%
  group_by(Site, Subject) %>%
  do(correlation_matrix = {
    df <- .[, -1] %>% select_if(is.numeric)
    if (ncol(df) > 1) {
      cor(df, use = "pairwise.complete.obs")
    } else {
      return(NA)
    }
  }) %>%
  ungroup() %>%
  slice(1)

cor_matrix_first_obs <- first_observation$correlation_matrix[[1]]

if (!all(is.na(cor_matrix_first_obs))) {
  melted_cor_matrix <- melt(cor_matrix_first_obs)
  print(ggplot(data = melted_cor_matrix, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1)) +
    coord_fixed() +
    labs(title = "Correlation Matrix for the First Observation",
         x = "Variables",
         y = "Variables"))
} else {
  print("No correlation matrix available for the first observation.")
}

##Find rows with N/A
na_rows <- apply(correlation_results_wide, 1, function(x) any(is.na(x)))
correlation_results_with_na <- correlation_results_wide[na_rows, ]
print(correlation_results_with_na)

})

## ================================================================================
## STAGE 4 — Merge Phenotype + Connectivity Data  (source: Merged_Data.R)
## Input: phenotype_data.csv + correlation_results_wide.csv  ->  Output: merged_data.csv
## ================================================================================
run_stage(4, "Merge Phenotype + Connectivity Data", {

n_data <- fread("phenotype_data.csv")
corr_data <- fread("correlation_results_wide.csv")

setDT(n_data)
setDT(corr_data)

dim(n_data) #[1] 568  15
dim(corr_data) #[1]  1053 76638

merged_data <- merge(
  n_data,
  corr_data,
  by.x = "SUB_ID",
  by.y = "Subject",
  all = FALSE #inner join
)

merged_data[, Subject := NULL]

dim(merged_data)#[1]   486 76652

write.csv(merged_data, "merged_data.csv")

})

## ================================================================================
## STAGE 5 — Principal Component Analysis  (source: PCA.R)
## Input: merged_data.csv  ->  Output: varimax_scores.csv, promax_scores.csv,
##        rotated loadings, scree plots, per-component loading CSVs
## ================================================================================
run_stage(5, "Principal Component Analysis", {

merged_data <- fread("merged_data.csv")

start_time <- Sys.time()

start_col <- which(names(merged_data) == "X.10_X.1")
end_col <- which(names(merged_data) == "X.99_X.98")

selected_vars <- merged_data[, start_col:end_col]
numeric_vars <- na.omit(as.data.frame(lapply(selected_vars, as.numeric)))

pca_result <- prcomp(numeric_vars, center = TRUE, scale. = TRUE)

eigenvalues <- pca_result$sdev^2
explained_var <- pca_result$sdev^2
prop_var <- explained_var / sum(explained_var)
cum_var <- cumsum(prop_var)

var_df <- data.frame(
  PC = paste0("PC", 1:length(prop_var)),
  Proportion = prop_var,
  Cumulative = cum_var,
  Explained_Variance = explained_var
)

print(var_df)
write.csv(var_df, "Explained Variance_PCA.csv")

explained_variance <- (eigenvalues / sum(eigenvalues)) * 100

png("scree_plot_explained_variance_25_intervals.png", width = 800, height = 600)
num_components <- length(explained_variance)
x_intervals <- seq(1, num_components, length.out = 25)
plot(explained_variance, type = "b",
     main = "Scree Plot - Explained Variance",
     xlab = "Principal Component",
     ylab = "Explained Variance (%)",
     pch = 19, col = "blue", xaxt = "n")
axis(1, at = x_intervals, labels = round(x_intervals))
abline(v = 20, col = "darkgreen", lty = 3, lwd = 2)
dev.off()
cat("scree_plot_explained_variance_25_intervals.png\n")

png("scree_plot_prcomp_25_intervals.png", width = 800, height = 600)
num_components <- length(eigenvalues)
x_intervals <- seq(1, num_components, length.out = 25)
plot(eigenvalues, type = "b", main = "Scree Plot", xlab = "Principal Component", ylab = "Eigenvalue",
     pch = 19, col = "blue", xaxt = "n")
axis(1, at = x_intervals, labels = round(x_intervals))
abline(h = 1, col = "red", lty = 2)
dev.off()
cat("Scree plot.png\n")

loadings <- as.matrix(pca_result$rotation)
loadings_50 <- as.matrix(pca_result$rotation[, 1:50])

varimax_result <- varimax(loadings_50)
rotated_loadings_varimax <- varimax_result$loadings

get_top_features <- function(rotated_loadings, top_n = 5) {
  rotated_loadings <- as.matrix(rotated_loadings)
  if (ncol(rotated_loadings) == 0 || nrow(rotated_loadings) == 0) {
    stop("Error: rotated_loadings has zero rows or columns.")
  }
  top_features_list <- list()
  for (i in seq_len(ncol(rotated_loadings))) {
    component <- rotated_loadings[, i]
    feature_names <- rownames(rotated_loadings)
    if (is.null(feature_names)) {
      stop("Error: Row names (features) are missing in rotated_loadings.")
    }
    component_features <- data.frame(
      Feature = feature_names,
      Loading = component
    )
    component_features <- component_features[order(abs(component_features$Loading), decreasing = TRUE), ]
    top_n <- min(top_n, nrow(component_features))
    top_features <- component_features[1:top_n, ]
    top_features_list[[paste0("Component_", i)]] <- top_features
  }
  return(top_features_list)
}

rotated_loadings_matrix <- as.matrix(rotated_loadings)
top_features_varimax <- get_top_features(rotated_loadings_matrix, top_n = 10)

# NOTE: relative subfolder (was an absolute path in the original script)
output_directory <- "Varimax Factor Importance"
if (!dir.exists(output_directory)) dir.create(output_directory)

for (i in seq_len(ncol(rotated_loadings))) {
  component_loadings <- rotated_loadings[, i, drop = FALSE]
  filename <- paste0(output_directory, "/component_", i, "_loadings.csv")
  write.csv(component_loadings, filename, row.names = TRUE)
  cat("Saved loadings for Component", i, "to", filename, "\n")
}

rotated_loadings <- as.matrix(rotated_loadings)
num_components <- min(20, ncol(rotated_loadings))
selected_loadings <- rotated_loadings[, 1:num_components]
loadings_dataframe <- as.data.frame(selected_loadings)
loadings_dataframe$Feature <- rownames(rotated_loadings)
loadings_dataframe <- loadings_dataframe[, c("Feature", colnames(loadings_dataframe)[-ncol(loadings_dataframe)])]
write.csv(loadings_dataframe, "rotated_loadings_20_components_varimax.csv", row.names = FALSE)

extreme_values <- data.frame(
  Feature = rownames(selected_loadings),
  Highest_Component = NA,
  Highest_Value = NA,
  Lowest_Component = NA,
  Lowest_Value = NA
)
for (i in 1:nrow(selected_loadings)) {
  variable_loadings <- selected_loadings[i, ]
  extreme_values$Highest_Component[i] <- which.max(variable_loadings)
  extreme_values$Highest_Value[i] <- max(variable_loadings)
  extreme_values$Lowest_Component[i] <- which.min(variable_loadings)
  extreme_values$Lowest_Value[i] <- min(variable_loadings)
}
write.csv(extreme_values, "Factor Importance_varimax.csv", row.names = FALSE)

numeric_vars_matrix <- as.matrix(numeric_vars)
varimax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings)

promax_result <- promax(loadings_50)
rotated_loadings_promax <- promax_result$loadings

rotated_loadings_matrix <- as.matrix(rotated_loadings_promax)
top_features_promax <- get_top_features(rotated_loadings_matrix, top_n = 10)

# NOTE: relative subfolder (was an absolute path in the original script)
output_directory <- "Promax Factor Importance"
if (!dir.exists(output_directory)) dir.create(output_directory)

for (i in seq_len(ncol(rotated_loadings_promax))) {
  component_loadings <- rotated_loadings_promax[, i, drop = FALSE]
  filename <- paste0(output_directory, "/component_", i, "_loadings.csv")
  write.csv(component_loadings, filename, row.names = TRUE)
  cat("Saved loadings for Component", i, "to", filename, "\n")
}

rotated_loadings_promax <- as.matrix(rotated_loadings_promax)
num_components <- min(20, ncol(rotated_loadings_promax))
selected_loadings <- rotated_loadings_promax[, 1:num_components]
loadings_dataframe <- as.data.frame(selected_loadings)
loadings_dataframe$Feature <- rownames(rotated_loadings_promax)
loadings_dataframe <- loadings_dataframe[, c("Feature", colnames(loadings_dataframe)[-ncol(loadings_dataframe)])]
write.csv(loadings_dataframe, "rotated_loadings_20_components_promax.csv", row.names = FALSE)

extreme_values <- data.frame(
  Feature = rownames(selected_loadings),
  Highest_Component = NA,
  Highest_Value = NA,
  Lowest_Component = NA,
  Lowest_Value = NA
)
for (i in 1:nrow(selected_loadings)) {
  variable_loadings <- selected_loadings[i, ]
  extreme_values$Highest_Component[i] <- which.max(variable_loadings)
  extreme_values$Highest_Value[i] <- max(variable_loadings)
  extreme_values$Lowest_Component[i] <- which.min(variable_loadings)
  extreme_values$Lowest_Value[i] <- min(variable_loadings)
}
write.csv(extreme_values, "Factor Importance_promax.csv", row.names = FALSE)

numeric_vars_matrix <- as.matrix(numeric_vars)
promax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings_promax)

rotated_loadings_varimax <- varimax_result$loadings
rotated_loadings_promax <- promax_result$loadings

threshold <- 0.00

get_variables_by_component_wide <- function(rotated_loadings, threshold) {
  results <- list()
  for (i in 1:ncol(rotated_loadings)) {
    component_vars <- rownames(rotated_loadings)[which(abs(rotated_loadings[, i]) > threshold)]
    results[[paste0("Component_", i)]] <- component_vars
  }
  max_length <- max(sapply(results, length))
  wide_data <- as.data.frame(do.call(cbind, lapply(results, function(x) {
    c(x, rep(NA, max_length - length(x)))
  })))
  colnames(wide_data) <- names(results)
  return(wide_data)
}

wide_varimax <- get_variables_by_component_wide(rotated_loadings_varimax, threshold)
wide_promax <- get_variables_by_component_wide(rotated_loadings_promax, threshold)

cat("\nVarimax Wide-Format Results:\n")
print(wide_varimax)
cat("\nPromax Wide-Format Results:\n")
print(wide_promax)

write.csv(varimax_scores, "varimax_scores.csv")
write.csv(promax_scores, "promax_scores.csv")
write.csv(wide_varimax, "varimax_pc_variables.csv")
write.csv(wide_promax, "promax_pc_variables.csv")

total_end_time <- Sys.time()
# FIX: original passed a string ("seconds\n") as round()'s digits arg, which
# errors at runtime. Split into two cat() calls instead.
cat("Total runtime:", round(difftime(total_end_time, start_time, units = "secs"), 2), "seconds\n")

})

## ================================================================================
## STAGE 6 — Dataset Creation  (source: Dataset_Creation.R)
## Input: merged_data.csv + varimax_scores.csv + promax_scores.csv
## Output: varimax_alldata.csv, promax_alldata.csv
## ================================================================================
run_stage(6, "Dataset Creation", {

merged_data <- fread("merged_data.csv")
varimax_scores <- fread("varimax_scores.csv")
promax_scores <- fread("promax_scores.csv")

alldata <- merged_data %>% select(SUB_ID, SITE_ID, DX_GROUP, DSM_IV_TR,
                                  AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM)

alldata$AVG_SOCIAL_INT <- round(alldata$AVG_SOCIAL * 10)
alldata$AVG_COMM_INT <- round(alldata$AVG_COMM * 10)

varimax_alldata <- combined_data <- cbind(alldata, varimax_scores)
varimax_alldata <- varimax_alldata %>% select(-V1)
promax_alldata <- combined_data <- cbind(alldata, promax_scores)
promax_alldata <- promax_alldata %>% select(-V1)

write.csv(varimax_alldata, "varimax_alldata.csv")
write.csv(promax_alldata, "promax_alldata.csv")

})

## ================================================================================
## STAGE 7 — Exploratory Data Analysis: PC Distributions  (source: Data Exploratory.R)
## Input: promax_alldata.csv, varimax_alldata.csv  ->  Output: KS_Test_Results.csv
## ================================================================================
run_stage(7, "EDA - PC Distributions (QQ / CDF / KS)", {

promax_data <- fread("promax_alldata.csv")
varimax_data <- fread("varimax_alldata.csv")

pc_columns <- paste0("PC", 1:20)

promax_data$DX_GROUP <- as.factor(promax_data$DX_GROUP)
varimax_data$DX_GROUP <- as.factor(varimax_data$DX_GROUP)

### QQ PLOT for PC1 to PC20 (Promax & Varimax) ###
generate_qq_plot <- function(data, dataset_name) {
  melted_data <- data %>%
    select(all_of(pc_columns)) %>%
    pivot_longer(cols = everything(), names_to = "Component", values_to = "Value")

  ggplot(melted_data, aes(sample = Value)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Component, scales = "free") +
    labs(title = paste("QQ Plot for", dataset_name, "Dataset"), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
}

qq_plot_promax <- generate_qq_plot(promax_data, "Promax")
qq_plot_varimax <- generate_qq_plot(varimax_data, "Varimax")

print(qq_plot_promax)
print(qq_plot_varimax)

### CDF PLOT for PC1 to PC20 by DX_GROUP ###
generate_cdf_plot <- function(data, dataset_name) {
  melted_data <- data %>%
    select(DX_GROUP, all_of(pc_columns)) %>%
    pivot_longer(cols = -DX_GROUP, names_to = "Component", values_to = "Value")

  ggplot(melted_data, aes(x = Value, color = DX_GROUP)) +
    stat_ecdf(size = 1.2) +
    facet_wrap(~ Component, scales = "free") +
    labs(title = paste("Cumulative Density Function (CDF) for", dataset_name, "Dataset"),
         x = "Principal Component Value", y = "Cumulative Probability",
         color = "DX_GROUP") +
    theme_minimal()
}

cdf_plot_promax <- generate_cdf_plot(promax_data, "Promax")
cdf_plot_varimax <- generate_cdf_plot(varimax_data, "Varimax")

print(cdf_plot_promax)
print(cdf_plot_varimax)

### KS TEST for PC1 to PC20 by DX_GROUP ###
perform_ks_test <- function(data, dataset_name) {
  ks_results <- data.frame(Variable = pc_columns, P_Value = NA)

  for (pc in pc_columns) {
    group1 <- data %>% filter(DX_GROUP == levels(DX_GROUP)[1]) %>% pull(pc)
    group2 <- data %>% filter(DX_GROUP == levels(DX_GROUP)[2]) %>% pull(pc)

    if (length(group1) > 0 & length(group2) > 0) {
      ks_test <- ks.test(group1, group2)
      ks_results$P_Value[ks_results$Variable == pc] <- ks_test$p.value
    }
  }

  ks_results$Dataset <- dataset_name
  return(ks_results)
}

ks_results_promax <- perform_ks_test(promax_data, "Promax")
ks_results_varimax <- perform_ks_test(varimax_data, "Varimax")

ks_results <- rbind(ks_results_promax, ks_results_varimax)
print(ks_results)
write.csv(ks_results, "KS_Test_Results.csv", row.names = FALSE)

})

## ================================================================================
## STAGE 8 — Exploratory Data Analysis: Demographics & PC Summaries  (source: EDA_initial.R)
## Input: merged_data.csv, varimax_alldata.csv, promax_alldata.csv
## ================================================================================
run_stage(8, "EDA - Demographics & PC Summaries", {

data <- fread("merged_data.csv")
varimax_alldata <- fread("varimax_alldata.csv")
promax_alldata <- fread("promax_alldata.csv")
data <- data %>% select(SITE_ID:AVG_COMM)

summary(data)

dx_group_table <- data %>%
  group_by(DX_GROUP) %>%
  summarize(
    Mean_Age = mean(AGE_AT_SCAN, na.rm = TRUE),
    SD_Age = sd(AGE_AT_SCAN, na.rm = TRUE),
    Male_Count = sum(SEX == 1, na.rm = TRUE),
    Female_Count = sum(SEX == 2, na.rm = TRUE)
  )
write.csv(dx_group_table, file = "dx_group_table.csv", row.names = FALSE)

age_plot <- ggplot(data, aes(x = AGE_AT_SCAN)) +
  geom_histogram(binwidth = 2, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of Age", x = "Age", y = "Frequency")
print(age_plot)
ggsave("age_plot.png", plot = age_plot)

dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP))) +
  geom_bar(fill = "green", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
  labs(title = "Distribution of Diagnostic Groups", x = "Diagnostic Group", y = "Count")
print(dx_plot)
ggsave("dx_plot.png", plot = dx_plot)

gender_plot <- ggplot(data, aes(x = as.factor(SEX))) +
  geom_bar(fill = "purple", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")
print(gender_plot)
ggsave("gender_plot.png", plot = gender_plot)

pairplot_vars <- data %>% select(DX_GROUP, AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM)
png("pairplot.png")
pairs(pairplot_vars, main = "Pairplot of Selected Key Variables")
dev.off()

age_dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP), y = AGE_AT_SCAN)) +
  geom_boxplot(fill = "orange", alpha = 0.7) +
  scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
  labs(title = "Age by Diagnostic Group", x = "Diagnostic Group", y = "Age")
print(age_dx_plot)
ggsave("age_dx_plot.png", plot = age_dx_plot)

social_age_plot <- ggplot(data, aes(x = AGE_AT_SCAN, y = AVG_SOCIAL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "AVG_SOCIAL vs Age", x = "Age", y = "AVG_SOCIAL")
print(social_age_plot)
ggsave("social_age_plot.png", plot = social_age_plot)

avg_social_sex_plot <- ggplot(data, aes(x = as.factor(SEX), y = AVG_SOCIAL)) +
  geom_boxplot(fill = "lightgreen", alpha = 0.7) +
  scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
  labs(title = "AVG_SOCIAL by Sex", x = "Sex", y = "AVG_SOCIAL")
print(avg_social_sex_plot)
ggsave("avg_social_sex_plot.png", plot = avg_social_sex_plot)

selected_vars <- c('ADI_R_SOCIAL_TOTAL_A', 'ADI_R_VERBAL_TOTAL_BV', 'ADOS_COMM',
                   'ADOS_SOCIAL', 'SRS_COGNITION', 'SRS_COMMUNICATION',
                   'AVG_SOCIAL', 'AVG_COMM')

summary_stats <- data %>%
  group_by(DX_GROUP) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(summary_stats)

for (var in selected_vars) {
  dx_plot <- ggplot(data, aes_string(x = "as.factor(DX_GROUP)", y = var)) +
    geom_boxplot(fill = "cyan", alpha = 0.7) +
    scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
    labs(title = paste(var, "by Diagnostic Group"), x = "Diagnostic Group", y = var)
  ggsave(paste0("dx_plot_", var, ".png"), plot = dx_plot)
}

anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ DX_GROUP")), data = data)
  anova_results[[var]] <- summary(model)
}
print(anova_results)

dx_group_pvalues <- sapply(anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(dx_group_pvalues), file = "dx_group_pvalues.csv", row.names = TRUE)

for (var in selected_vars) {
  trend_plot <- ggplot(data, aes_string(x = "AGE_AT_SCAN", y = var, color = "as.factor(DX_GROUP)")) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste(var, "Trends by Age and Diagnostic Group"), x = "Age", y = var, color = "Diagnostic Group")
  print(trend_plot)
}

sex_summary_stats <- data %>%
  group_by(SEX) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(sex_summary_stats)

for (var in selected_vars) {
  sex_plot <- ggplot(data, aes_string(x = "as.factor(SEX)", y = var)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
    labs(title = paste(var, "by Sex"), x = "Sex", y = var)
  ggsave(paste0("sex_plot_", var, ".png"), plot = sex_plot)
}

sex_anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ SEX")), data = data)
  sex_anova_results[[var]] <- summary(model)
}
print(sex_anova_results)
sex_anova_pvalues <- sapply(sex_anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(sex_anova_pvalues), file = "sex_anova_pvalues.csv", row.names = TRUE)

data <- data %>% mutate(AGE_INTERVAL = cut(AGE_AT_SCAN,
                                           breaks = c(-Inf, 10, 20, 30, 40, Inf),
                                           labels = c("<10", "10-20", "20-30", "30-40", ">40")))

AGE_INTERVAL_summary_stats <- data %>%
  group_by(AGE_INTERVAL) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(AGE_INTERVAL_summary_stats)

for (var in selected_vars) {
  AGE_INTERVAL_plot <- ggplot(data, aes(x = AGE_INTERVAL, y = .data[[var]])) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    labs(title = paste(var, "by AGE_INTERVAL"), x = "AGE_INTERVAL", y = var)
  ggsave(paste0("AGE_INTERVAL_plot_", var, ".png"), plot = AGE_INTERVAL_plot)
}

AGE_INTERVAL_anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ AGE_INTERVAL")), data = data)
  AGE_INTERVAL_anova_results[[var]] <- summary(model)
}
print(AGE_INTERVAL_anova_results)
age_interval_pvalues <- sapply(AGE_INTERVAL_anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(age_interval_pvalues), file = "age_interval_anova_pvalues.csv", row.names = TRUE)

age_int_grouped_table <- data %>%
  group_by(AGE_INTERVAL) %>%
  summarize(across(ADI_R_SOCIAL_TOTAL_A:AVG_COMM, mean, na.rm = TRUE))
write.csv(age_int_grouped_table, file = "grouped_table_by_age_interval.csv", row.names = FALSE)

varimax_data <- varimax_alldata %>% select(where(is.numeric)) %>% select(PC1:PC20)
promax_data <- promax_alldata %>% select(where(is.numeric)) %>% select(PC1:PC20)

varimax_corr <- cor(varimax_data, use = "complete.obs")
promax_corr <- cor(promax_data, use = "complete.obs")

png("varimax_correlation_matrix.png", width = 800, height = 800)
corrplot(varimax_corr, method = "color", title = "Correlation Matrix: Varimax", tl.cex = 0.8, cl.cex = 0.8, addgrid.col = NA)
dev.off()

png("promax_correlation_matrix.png", width = 800, height = 800)
corrplot(promax_corr, method = "color", title = "Correlation Matrix: Promax", tl.cex = 0.8, cl.cex = 0.8, addgrid.col = NA)
dev.off()

##PC Average for DX_GROUP
promax_alldata <- read.csv("promax_alldata.csv")
varimax_alldata <- read.csv("varimax_alldata.csv")

summarize_pcs_by_group <- function(data, dataset_name) {
  data %>%
    select(DX_GROUP, PC1:PC20) %>%
    group_by(DX_GROUP) %>%
    summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-DX_GROUP, names_to = "PC", values_to = "Mean") %>%
    pivot_wider(names_from = DX_GROUP, values_from = Mean) %>%
    rename_with(~ paste0(dataset_name, "_", .), -PC)
}

promax_summary <- summarize_pcs_by_group(promax_alldata, "Promax")
varimax_summary <- summarize_pcs_by_group(varimax_alldata, "Varimax")

combined_summary <- full_join(promax_summary, varimax_summary, by = "PC")
write.csv(combined_summary, file = "pc_summary_by_dx_group.csv", row.names = FALSE)

visualization_data <- combined_summary %>%
  pivot_longer(cols = -PC, names_to = "Dataset_DX_GROUP", values_to = "Mean") %>%
  separate(Dataset_DX_GROUP, into = c("Dataset", "DX_GROUP"), sep = "_") %>%
  mutate(DX_GROUP = ifelse(DX_GROUP == "1", "Autism", "Control")) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:20)))

pc_summary_plot <- ggplot(visualization_data, aes(x = PC, y = Mean, fill = DX_GROUP)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Dataset, ncol = 1) +
  theme_minimal() +
  labs(
    title = "PC Summary Statistics by DX_GROUP",
    x = "Principal Component (PC)",
    y = "Mean Value",
    fill = "DX_GROUP"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("pc_summary_plot.png", plot = pc_summary_plot, width = 10, height = 8)

prepare_cdf_data <- function(data, dataset_name) {
  data %>%
    select(starts_with("PC")) %>%
    pivot_longer(cols = everything(), names_to = "PC", values_to = "Value") %>%
    filter(!is.na(Value) & PC %in% paste0("PC", 1:20)) %>%
    group_by(PC) %>%
    mutate(CDF = ecdf(Value)(Value)) %>%
    ungroup() %>%
    mutate(Dataset = dataset_name)
}

promax_cdf_data <- prepare_cdf_data(promax_alldata, "Promax")
varimax_cdf_data <- prepare_cdf_data(varimax_alldata, "Varimax")

cdf_data <- bind_rows(promax_cdf_data, varimax_cdf_data) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:20)))

ks_test_results <- lapply(paste0("PC", 1:20), function(pc) {
  promax_values <- promax_alldata[[pc]]
  varimax_values <- varimax_alldata[[pc]]
  ks_result <- ks.test(promax_values, varimax_values, alternative = "two.sided")
  data.frame(
    PC = pc,
    D_statistic = ks_result$statistic,
    P_value = ks_result$p.value
  )
})

ks_test_summary <- do.call(rbind, ks_test_results)
write.csv(ks_test_summary, file = "ks_test_results.csv", row.names = FALSE)
print(ks_test_summary)

cdf_plot <- ggplot(cdf_data, aes(x = Value, y = CDF, color = Dataset)) +
  geom_line(size = 1.2, alpha = 0.8) +
  facet_wrap(~ PC, ncol = 4, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Cumulative Distribution Functions of Principal Components",
    x = "Value",
    y = "Cumulative Probability",
    color = "Dataset"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave("cdf_plot_pc_variables.png", plot = cdf_plot, width = 12, height = 10)
print(cdf_plot)

})

## ================================================================================
## STAGE 9 — PC4 Brain Network Visualization: Varimax  (source: PC4_Varimax_Visual.R)
## Input: Varimax Factor Importance/component_4_loadings.csv, ROI_LABELS_FILE
## ================================================================================
run_stage(9, "PC4 Brain Network Visualization (Varimax)", {

loadings <- read_csv(file.path("Varimax Factor Importance", "component_4_loadings.csv"), col_names = FALSE)
roi_labels <- read_csv(ROI_LABELS_FILE, col_names = FALSE)

colnames(loadings) <- c("ROI_Connection", "Loading_Value")
loadings$Loading_Value <- as.numeric(loadings$Loading_Value)

loadings <- loadings %>%
  mutate(
    ROI_1 = sub("_.*", "", ROI_Connection),
    ROI_2 = sub(".*_", "", ROI_Connection)
  )

loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)

roi_labels <- roi_labels %>%
  select(ROI_Number = 1, AAL_Label = 5)

loadings <- loadings %>%
  left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
  rename(ROI_1_Label = AAL_Label) %>%
  left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
  rename(ROI_2_Label = AAL_Label)

loadings <- loadings %>%
  mutate(Connection_Label = paste(ROI_1_Label, "to", ROI_2_Label))

top20 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:20)

edges <- top20 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

nodes <- data.frame(name = unique(c(edges$from, edges$to)))

nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -1.5, max = -0.5), runif(n()/2, min = 0.5, max = 1.5)),
    y = seq(-1, 1, length.out = n())
  )

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>% mutate(x = nodes$x, y = nodes$y)

print(ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8) +
  geom_node_point(size = 6, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = "Top 20 ROI Connections (Brain Hemisphere Layout)- Varimax",
       edge_width = "Strength",
       edge_color = "Loading Value") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

top20_rois <- unique(c(top20$ROI_1_Label, top20$ROI_2_Label))
write.csv(top20, "Top20_Varimax.csv")

roi_brain_region <- function(label) {
  if (grepl("Frontal", label)) return("Frontal Lobe")
  else if (grepl("Temporal", label)) return("Temporal Lobe")
  else if (grepl("Parietal", label)) return("Parietal Lobe")
  else if (grepl("Occipital", label)) return("Occipital Lobe")
  else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
  else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
  else return("Other")
}

top20_long <- top20 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

brain_heat_data <- top20_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

})

## ================================================================================
## STAGE 10 — PC4 Brain Network Visualization: Promax  (source: PC4_Promax_Visual.R)
## FLAG: original reads "component_4_loadings_promax.csv", but Stage 5 (PCA.R)
## writes "component_4_loadings.csv" into the Promax Factor Importance folder —
## the two names never matched in the original scripts either. Reading the file
## that actually exists (component_4_loadings.csv) below; rename/adjust if you
## intended a differently-named file.
## ================================================================================
run_stage(10, "PC4 Brain Network Visualization (Promax)", {

loadings <- read_csv(file.path("Promax Factor Importance", "component_4_loadings.csv"), col_names = FALSE)
roi_labels <- read_csv(ROI_LABELS_FILE, col_names = FALSE)

colnames(loadings) <- c("ROI_Connection", "Loading_Value")
loadings$Loading_Value <- as.numeric(loadings$Loading_Value)

loadings <- loadings %>%
  mutate(
    ROI_1 = sub("_.*", "", ROI_Connection),
    ROI_2 = sub(".*_", "", ROI_Connection)
  )

loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)

roi_labels <- roi_labels %>%
  select(ROI_Number = 1, AAL_Label = 5)

loadings <- loadings %>%
  left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
  rename(ROI_1_Label = AAL_Label) %>%
  left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
  rename(ROI_2_Label = AAL_Label)

loadings <- loadings %>%
  mutate(Connection_Label = paste(ROI_1_Label, "to", ROI_2_Label))

top20 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:20)

edges <- top20 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

nodes <- data.frame(name = unique(c(edges$from, edges$to)))

nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -1.5, max = -0.5), runif(n()/2, min = 0.5, max = 1.5)),
    y = seq(-1, 1, length.out = n())
  )

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>% mutate(x = nodes$x, y = nodes$y)

print(ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8) +
  geom_node_point(size = 6, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = "Top 20 ROI Connections (Brain Hemisphere Layout)- Promax",
       edge_width = "Strength",
       edge_color = "Loading Value") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

top20_rois <- unique(c(top20$ROI_1_Label, top20$ROI_2_Label))
write.csv(top20, "Top20_Promax.csv")

roi_brain_region <- function(label) {
  if (grepl("Frontal", label)) return("Frontal Lobe")
  else if (grepl("Temporal", label)) return("Temporal Lobe")
  else if (grepl("Parietal", label)) return("Parietal Lobe")
  else if (grepl("Occipital", label)) return("Occipital Lobe")
  else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
  else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
  else return("Other")
}

top20_long <- top20 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

brain_heat_data <- top20_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

})

## ================================================================================
## STAGE 11 — PC17 Communication Brain Network Visualization  (source: PC4_Promax_Visual_comm.R)
## FLAG: same filename-mismatch issue as Stage 10 — original reads
## "component_17_loadings_comm.csv" which PCA.R never writes under that name;
## reading "component_17_loadings.csv" from the Promax folder instead.
## ================================================================================
run_stage(11, "PC17 Communication Brain Network Visualization", {

loadings <- read_csv(file.path("Promax Factor Importance", "component_17_loadings.csv"), col_names = FALSE)
roi_labels <- read_csv(ROI_LABELS_FILE, col_names = FALSE)

colnames(loadings) <- c("ROI_Connection", "Loading_Value")
loadings$Loading_Value <- as.numeric(loadings$Loading_Value)

loadings <- loadings %>%
  mutate(
    ROI_1 = sub("_.*", "", ROI_Connection),
    ROI_2 = sub(".*_", "", ROI_Connection)
  )

loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)

roi_labels <- roi_labels %>%
  select(ROI_Number = 1, AAL_Label = 5)

loadings <- loadings %>%
  left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
  rename(ROI_1_Label = AAL_Label) %>%
  left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
  rename(ROI_2_Label = AAL_Label)

loadings <- loadings %>%
  mutate(Connection_Label = paste(ROI_1_Label, "to", ROI_2_Label))

top30 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:30)

write.csv(top30, "Top30_Comm.csv")

edges <- top30 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

nodes <- data.frame(name = unique(c(edges$from, edges$to)))

nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -2, max = -0.8), runif(n()/2, min = 0.8, max = 2)),
    y = seq(-2, 2, length.out = n())
  )

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>% mutate(x = nodes$x, y = nodes$y)

print(ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.6) +
  geom_node_point(size = 4, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "plain") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.3, 1.5)) +
  theme_void() +
  labs(
    title = "Top 30 ROI Connections (Brain Hemisphere Layout) - Communication (PC17)",
    edge_width = "Strength",
    edge_color = "Loading Value"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

top30_rois <- unique(c(top30$ROI_1_Label, top30$ROI_2_Label))
write.csv(top30, "top30_Promax.csv")

roi_brain_region <- function(label) {
  if (grepl("Frontal", label)) return("Frontal Lobe")
  else if (grepl("Temporal", label)) return("Temporal Lobe")
  else if (grepl("Parietal", label)) return("Parietal Lobe")
  else if (grepl("Occipital", label)) return("Occipital Lobe")
  else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
  else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
  else return("Other")
}

top30_long <- top30 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

brain_heat_data <- top30_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

#Get Median Data of Top30 FC by DX_GROUP
merged_data <- fread("merged_data.csv")

roi_vars <- top30$ROI_Connection
roi_vars <- roi_vars[roi_vars %in% names(merged_data)]

if (!"DX_GROUP" %in% names(merged_data)) {
  stop("DX_GROUP column not found in merged_data.")
}

roi_data <- merged_data[, c("DX_GROUP", roi_vars), with = FALSE]

median_by_group <- roi_data %>%
  group_by(DX_GROUP) %>%
  summarise(across(all_of(roi_vars), median, na.rm = TRUE), .groups = "drop")

print(median_by_group)
write.csv(median_by_group, "median_roi_by_dx_group_Comm.csv", row.names = FALSE)

median_long <- median_by_group %>%
  pivot_longer(cols = -DX_GROUP, names_to = "ROI_Connection", values_to = "Median_Value")

print(ggplot(median_long, aes(x = ROI_Connection, y = Median_Value, fill = as.factor(DX_GROUP))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(
    values = c("1" = "#1f77b4", "2" = "#ff7f0e"),
    name = "DX Group",
    labels = c("1" = "Autism", "2" = "Control")
  ) +
  labs(
    title = "Median ROI Connectivity by Diagnostic Group",
    x = "ROI Connection",
    y = "Median Value"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  ))

})

## ================================================================================
## STAGE 12 — PC15 Social Cognition Brain Network Visualization  (source: PC4_Promax_Visual_social.R)
## FLAG: same filename-mismatch issue as Stage 10/11 — reading
## "component_15_loadings.csv" from the Promax folder instead of the original
## "component_15_loadings_soc.csv" (which PCA.R never produced under that name).
## ================================================================================
run_stage(12, "PC15 Social Cognition Brain Network Visualization", {

loadings <- read_csv(file.path("Promax Factor Importance", "component_15_loadings.csv"), col_names = FALSE)
roi_labels <- read_csv(ROI_LABELS_FILE, col_names = FALSE)

colnames(loadings) <- c("ROI_Connection", "Loading_Value")
loadings$Loading_Value <- as.numeric(loadings$Loading_Value)

loadings <- loadings %>%
  mutate(
    ROI_1 = sub("_.*", "", ROI_Connection),
    ROI_2 = sub(".*_", "", ROI_Connection)
  )

loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)

roi_labels <- roi_labels %>%
  select(ROI_Number = 1, AAL_Label = 5)

loadings <- loadings %>%
  left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
  rename(ROI_1_Label = AAL_Label) %>%
  left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
  rename(ROI_2_Label = AAL_Label)

loadings <- loadings %>%
  mutate(Connection_Label = paste(ROI_1_Label, "to", ROI_2_Label))

top30 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:30)

edges <- top30 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

nodes <- data.frame(name = unique(c(edges$from, edges$to)))

nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -2, max = -0.8), runif(n()/2, min = 0.8, max = 2)),
    y = seq(-2, 2, length.out = n())
  )

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>% mutate(x = nodes$x, y = nodes$y)

print(ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.6) +
  geom_node_point(size = 4, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "plain") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.3, 1.5)) +
  theme_void() +
  labs(
    title = "Top 30 ROI Connections (Brain Hemisphere Layout) - Social Cognition (PC15)",
    edge_width = "Strength",
    edge_color = "Loading Value"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

top30_rois <- unique(c(top30$ROI_1_Label, top30$ROI_2_Label))
write.csv(top30, "top30_social.csv")

roi_brain_region <- function(label) {
  if (grepl("Frontal", label)) return("Frontal Lobe")
  else if (grepl("Temporal", label)) return("Temporal Lobe")
  else if (grepl("Parietal", label)) return("Parietal Lobe")
  else if (grepl("Occipital", label)) return("Occipital Lobe")
  else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
  else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
  else return("Other")
}

top30_long <- top30 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

brain_heat_data <- top30_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

#Get Median Data of Top30 FC by DX_GROUP
merged_data <- fread("merged_data.csv")

roi_vars <- top30$ROI_Connection
roi_vars <- roi_vars[roi_vars %in% names(merged_data)]

if (!"DX_GROUP" %in% names(merged_data)) {
  stop("DX_GROUP column not found in merged_data.")
}

roi_data <- merged_data[, c("DX_GROUP", roi_vars), with = FALSE]

median_by_group <- roi_data %>%
  group_by(DX_GROUP) %>%
  summarise(across(all_of(roi_vars), median, na.rm = TRUE), .groups = "drop")

print(median_by_group)
write.csv(median_by_group, "median_roi_by_dx_group_SOC.csv", row.names = FALSE)

median_long <- median_by_group %>%
  pivot_longer(cols = -DX_GROUP, names_to = "ROI_Connection", values_to = "Median_Value")

print(ggplot(median_long, aes(x = ROI_Connection, y = Median_Value, fill = as.factor(DX_GROUP))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(
    values = c("1" = "#1f77b4", "2" = "#ff7f0e"),
    name = "DX Group",
    labels = c("1" = "Autism", "2" = "Control")
  ) +
  labs(
    title = "Median ROI Connectivity by Diagnostic Group",
    x = "ROI Connection",
    y = "Median Value"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  ))

})

## ================================================================================
## STAGE 13 — MANCOVA Assumption Testing: Communication  (source: MANCOVA_Assumption Test_comm.R)
## ================================================================================
run_stage(13, "MANCOVA Assumption Testing - Communication", {

pdata <- read.csv("promax_alldata.csv")
vdata <- read.csv("varimax_alldata.csv")

dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, dv_names]
dependent_vars_v <- vdata[, dv_names]

pdata$DX_GROUP <- as.factor(pdata$DX_GROUP)
vdata$DX_GROUP <- as.factor(vdata$DX_GROUP)
independent_var <- "DX_GROUP"

covariates <- c("AVG_COMM_INT", "AGE_AT_SCAN", "SEX")

vif_test_p <- vif(lm(PC1 ~ DX_GROUP + AVG_COMM_INT + AGE_AT_SCAN + SEX, data = pdata))
print(vif_test_p)
vif_test_v <- vif(lm(PC1 ~ DX_GROUP + AVG_COMM_INT + AGE_AT_SCAN + SEX, data = vdata))
print(vif_test_v)

mardia_result <- mvn(dependent_vars_p, mvnTest = "mardia")
mardia_result_v <- mvn(dependent_vars_v, mvnTest = "mardia")
print(mardia_result$multivariateNormality)
print(mardia_result_v$multivariateNormality)

levene_results <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = pdata)
})
names(levene_results) <- dv_names
print(levene_results)

levene_results_v <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = vdata)
})
names(levene_results_v) <- dv_names
print(levene_results_v)

box_test <- boxM(dependent_vars_p, pdata[[independent_var]])
box_test_v <- boxM(dependent_vars_v, vdata[[independent_var]])
print(box_test)
print(box_test_v)

plot_data <- pdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list <- lapply(covariates, function(cov) {
  ggplot(plot_data, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables (Promax)"))
})
print(plot_list)

plot_data_v <- vdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list_v <- lapply(covariates, function(cov) {
  ggplot(plot_data_v, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables (Varimax)"))
})
print(plot_list_v)

mancova_model <- manova(as.matrix(dependent_vars_p) ~ DX_GROUP * AVG_COMM_INT +
                          AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                          SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                        data = pdata)

mancova_model_v <- manova(as.matrix(dependent_vars_v) ~ DX_GROUP * AVG_COMM_INT +
                            AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                            SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                          data = vdata)

anova_mancova <- summary.aov(mancova_model)
anova_mancova_v <- summary.aov(mancova_model_v)
print(anova_mancova)
print(anova_mancova_v)

})

## ================================================================================
## STAGE 14 — MANCOVA Assumption Testing: Social  (source: MANCOVA_Assumption Test_social.R)
## ================================================================================
run_stage(14, "MANCOVA Assumption Testing - Social", {

pdata <- read.csv("promax_alldata.csv")
vdata <- read.csv("varimax_alldata.csv")

dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, dv_names]
dependent_vars_v <- vdata[, dv_names]

pdata$DX_GROUP <- as.factor(pdata$DX_GROUP)
vdata$DX_GROUP <- as.factor(vdata$DX_GROUP)
independent_var <- "DX_GROUP"

covariates <- c("AVG_SOCIAL_INT", "AGE_AT_SCAN", "SEX")

vif_test_p <- vif(lm(PC1 ~ DX_GROUP + AVG_SOCIAL_INT + AGE_AT_SCAN + SEX, data = pdata))
print(vif_test_p)
vif_test_v <- vif(lm(PC1 ~ DX_GROUP + AVG_SOCIAL_INT + AGE_AT_SCAN + SEX, data = vdata))
print(vif_test_v)

mardia_result <- mvn(dependent_vars_p, mvnTest = "mardia")
mardia_result_v <- mvn(dependent_vars_v, mvnTest = "mardia")
print(mardia_result$multivariateNormality)
print(mardia_result_v$multivariateNormality)

levene_results <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = pdata)
})
names(levene_results) <- dv_names
print(levene_results)

levene_results_v <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = vdata)
})
names(levene_results_v) <- dv_names
print(levene_results_v)

box_test <- boxM(dependent_vars_p, pdata[[independent_var]])
box_test_v <- boxM(dependent_vars_v, vdata[[independent_var]])
print(box_test)
print(box_test_v)

plot_data <- pdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list <- lapply(covariates, function(cov) {
  ggplot(plot_data, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables"))
})
print(plot_list)

plot_data_v <- vdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list_v <- lapply(covariates, function(cov) {
  ggplot(plot_data_v, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables (Varimax)"))
})
print(plot_list_v)

mancova_model <- manova(as.matrix(dependent_vars_p) ~ DX_GROUP * AVG_SOCIAL_INT +
                          AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT +
                          SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                        data = pdata)

mancova_model_v <- manova(as.matrix(dependent_vars_v) ~ DX_GROUP * AVG_SOCIAL_INT +
                            AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT +
                            SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                          data = vdata)

anova_mancova <- summary.aov(mancova_model)
anova_mancova_v <- summary.aov(mancova_model_v)
print(anova_mancova)
print(anova_mancova_v)

})

## ================================================================================
## STAGE 15 — MANCOVA: Communication  (source: MANCOVA_alldata_comm.R)
## ================================================================================
run_stage(15, "MANCOVA - Communication", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_COMM_INT +
        SEX*DX_GROUP + SEX*AVG_COMM_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "varimax_detailed_summary_df[comm].csv")
write.csv(varimax_summary_df, "varimax_summary_df[comm].csv")
write.csv(promax_detailed_summary_df, "promax_detailed_summary_df[comm].csv")
write.csv(promax_summary_df, "promax_summary_df[comm].csv")

})

## ================================================================================
## STAGE 16 — MANCOVA: Social  (source: MANCOVA_alldata_social.R)
## ================================================================================
run_stage(16, "MANCOVA - Social", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_SOCIAL_INT +
        SEX*DX_GROUP + SEX*AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "varimax_detailed_summary_df[social].csv")
write.csv(varimax_summary_df, "varimax_summary_df[social].csv")
write.csv(promax_detailed_summary_df, "promax_detailed_summary_df[social].csv")
write.csv(promax_summary_df, "promax_summary_df[social].csv")

})

## ================================================================================
## STAGE 17 — Permutation MANCOVA: Communication  (source: Perm_MANCOVA_alldata_comm.R)
## ================================================================================
run_stage(17, "Permutation MANCOVA - Communication", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

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

})

## ================================================================================
## STAGE 18 — Permutation MANCOVA: Social  (source: Perm_MANCOVA_alldata_social.R)
## ================================================================================
run_stage(18, "Permutation MANCOVA - Social", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, ..dv_names]
dependent_vars_v <- vdata[, ..dv_names]

set.seed(123)
perm_mancova_p <- anova.cca(capscale(as.matrix(dependent_vars_p) ~ DX_GROUP * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = pdata), by = "margin", permutations = 10000)
print(perm_mancova_p)
write.csv(perm_mancova_p, "PermMANCOVA_alldata_promax_social.csv")

perm_mancova_v <- anova.cca(capscale(as.matrix(dependent_vars_v) ~ DX_GROUP * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = vdata), by = "margin", permutations = 10000)
print(perm_mancova_v)
write.csv(perm_mancova_v, "PermMANCOVA_alldata_varimax_social.csv")

})

## ================================================================================
## STAGE 19 — SVM VIF Multicollinearity Test  (source: SVM_VIF_Test.R)
## ================================================================================
run_stage(19, "SVM VIF Multicollinearity Test", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

predictor_vars <- paste0("PC", 1:20)

vdata[, (predictor_vars) := lapply(.SD, as.numeric), .SDcols = predictor_vars]
pdata[, (predictor_vars) := lapply(.SD, as.numeric), .SDcols = predictor_vars]

lm_vdata <- lm(DX_GROUP ~ ., data = vdata[, c("DX_GROUP", ..predictor_vars), with = FALSE])
lm_pdata <- lm(DX_GROUP ~ ., data = pdata[, c("DX_GROUP", ..predictor_vars), with = FALSE])

vif_vdata <- vif(lm_vdata)
vif_pdata <- vif(lm_pdata)

cat("VIF Results for Varimax Dataset:\n")
print(vif_vdata)
cat("\nVIF Results for Promax Dataset:\n")
print(vif_pdata)

write.csv(data.frame(Variable = names(vif_vdata), VIF = vif_vdata), "vif_results_varimax.csv", row.names = FALSE)
write.csv(data.frame(Variable = names(vif_pdata), VIF = vif_pdata), "vif_results_promax.csv", row.names = FALSE)

})

## ================================================================================
## STAGE 20 — SVM Classification Model  (source: SVM.R)
## Output: testing/training splits, svm_oos_output_*.csv, svm_oos_output_*_all.csv,
##         confusion matrices, ROC curves
## FIX: original had a stray trailing comma in the vdata out-of-sample roc() call
## ("svm_vdata_predictions[, 2],  )") which is an R syntax error and would have
## prevented this whole file — and therefore the whole combined script — from
## parsing. Removed the trailing comma below.
## ================================================================================
run_stage(20, "SVM Classification Model", {

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

##varimax dataset training and test
vdata$DX_GROUP <- as.factor(vdata$DX_GROUP)
set.seed(123)

vdata <- vdata %>% mutate(RowID = row_number())

testing_vdata <- vdata %>%
  filter(DX_GROUP %in% c("1", "2")) %>%
  group_by(DX_GROUP) %>%
  slice_sample(prop = 0.1) %>%
  ungroup()

training_vdata <- vdata %>%
  filter(!RowID %in% testing_vdata$RowID)

testing_vdata <- testing_vdata %>% select(-RowID)
training_vdata <- training_vdata %>% select(-RowID)

write.csv(testing_vdata, "testing_vdata.csv")
write.csv(training_vdata, "training_vdata.csv")

##promax dataset training and test
pdata$DX_GROUP <- as.factor(pdata$DX_GROUP)
set.seed(123)

pdata <- pdata %>% mutate(RowID = row_number())

testing_pdata <- pdata %>%
  filter(DX_GROUP %in% c("1", "2")) %>%
  group_by(DX_GROUP) %>%
  slice_sample(prop = 0.1) %>%
  ungroup()

training_pdata <- pdata %>%
  filter(!RowID %in% testing_pdata$RowID)

testing_pdata <- testing_pdata %>% select(-RowID)
training_pdata <- training_pdata %>% select(-RowID)

write.csv(testing_pdata, "testing_pdata.csv")
write.csv(training_pdata, "training_pdata.csv")

##SVM implementation
training_vdata$DX_GROUP <- as.factor(training_vdata$DX_GROUP)
training_pdata$DX_GROUP <- as.factor(training_pdata$DX_GROUP)

levels(training_vdata$DX_GROUP) <- make.names(levels(training_vdata$DX_GROUP))
levels(training_pdata$DX_GROUP) <- make.names(levels(training_pdata$DX_GROUP))

predictors <- paste0("PC", 1:20)
dependent_var <- "DX_GROUP"

custom_auc <- function(data, lev = NULL, model = NULL) {
  roc_obj <- roc(data$obs, as.numeric(data[, "pred"]))
  auc_value <- auc(roc_obj)
  return(c(AUC = as.numeric(auc_value)))
}

train_svm <- function(training_data) {
  training_data$DX_GROUP <- as.factor(training_data$DX_GROUP)

  train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    savePredictions = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = TRUE
  )

  svm_grid <- expand.grid(
    C = 2^(-5:5),
    sigma = 2^(-5:5)
  )

  svm_model <- train(
    DX_GROUP ~ .,
    data = training_data[, c("DX_GROUP", paste0("PC", 1:20))],
    method = "svmRadial",
    metric = "ROC",
    probability = TRUE,
    trControl = train_control,
    tuneGrid = svm_grid
  )

  return(svm_model)
}

cat("Training SVM for training_vdata...\n")
svm_vdata <- train_svm(training_vdata)

cat("\nTraining SVM for training_pdata...\n")
svm_pdata <- train_svm(training_pdata)

cat("\nSVM Model Summary for training_vdata:\n")
print(svm_vdata) #sigma = 0.03125 and C = 8

cat("\nSVM Model Summary for training_pdata:\n")
print(svm_pdata) #sigma = 0.03125 and C = 8

roc_vdata <- roc(
  response = svm_vdata$pred$obs,
  predictor = svm_vdata$pred$X1,
  levels = rev(levels(svm_vdata$pred$obs))
)
cat("\nAUC for training_vdata SVM model:\n")
print(auc(roc_vdata)) #Area under the curve: 0.6915 CV - Insample

roc_pdata <- roc(
  response = svm_pdata$pred$obs,
  predictor = svm_pdata$pred$X1,
  levels = rev(levels(svm_pdata$pred$obs))
)
cat("\nAUC for training_pdata SVM model:\n")
print(auc(roc_pdata)) #Area under the curve: 0.6997 CV - Insample

cat("Applying SVM to training_vdata...\n")
svm_vdata_predictions_all <- predict(svm_vdata, newdata = vdata, type = "prob")
svm_vdata_classes_all <- predict(svm_vdata, newdata = vdata, type = "raw")

svm_vdata_classes_all <- as.character(svm_vdata_classes_all)
svm_vdata_classes_all[svm_vdata_classes_all == "X1"] <- "1"
svm_vdata_classes_all[svm_vdata_classes_all == "X2"] <- "2"
svm_vdata_classes_all <- factor(svm_vdata_classes_all, levels = levels(vdata$DX_GROUP))

cat("Applying SVM to training_pdata...\n")
svm_pdata_predictions_all <- predict(svm_pdata, newdata = pdata, type = "prob")
svm_pdata_classes_all <- predict(svm_pdata, newdata = pdata, type = "raw")

svm_pdata_classes_all <- as.character(svm_pdata_classes_all)
svm_pdata_classes_all[svm_pdata_classes_all == "X1"] <- "1"
svm_pdata_classes_all[svm_pdata_classes_all == "X2"] <- "2"
svm_pdata_classes_all <- factor(svm_pdata_classes_all, levels = levels(pdata$DX_GROUP))

predictions_vdata_all <- as.data.frame(svm_vdata_classes_all)
predictions_pdata_all <- as.data.frame(svm_pdata_classes_all)
svm_oos_output_vdata_all <- cbind(vdata, predictions_vdata_all)
svm_oos_output_pdata_all <- cbind(pdata, predictions_pdata_all)

write.csv(svm_oos_output_vdata_all, "svm_oos_output_vdata_all.csv")
write.csv(svm_oos_output_pdata_all, "svm_oos_output_pdata_all.csv")

##Performance evaluation in OOS data
cat("Applying SVM to testing_vdata...\n")
svm_vdata_predictions <- predict(svm_vdata, newdata = testing_vdata, type = "prob")
svm_vdata_classes <- predict(svm_vdata, newdata = testing_vdata, type = "raw")

roc_vdata <- roc(
  testing_vdata$DX_GROUP,
  svm_vdata_predictions[, 2]
)
auc_vdata <- auc(roc_vdata)
cat("\nAUC for testing_vdata SVM model:\n")
print(auc_vdata) #Area under the curve: 0.8056

svm_vdata_classes <- as.character(svm_vdata_classes)
svm_vdata_classes[svm_vdata_classes == "X1"] <- "1"
svm_vdata_classes[svm_vdata_classes == "X2"] <- "2"
svm_vdata_classes <- factor(svm_vdata_classes, levels = levels(testing_vdata$DX_GROUP))

testing_vdata$DX_GROUP <- factor(testing_vdata$DX_GROUP, levels = c("1", "2"))
svm_vdata_classes <- factor(svm_vdata_classes, levels = levels(testing_vdata$DX_GROUP))
conf_matrix_vdata <- confusionMatrix(svm_vdata_classes, testing_vdata$DX_GROUP)
cat("\nConfusion Matrix for testing_vdata SVM model:\n")
print(conf_matrix_vdata)

cat("\nApplying SVM to testing_pdata...\n")
svm_pdata_predictions <- predict(svm_pdata, newdata = testing_pdata, type = "prob")
svm_pdata_classes <- predict(svm_pdata, newdata = testing_pdata, type = "raw")

svm_pdata_classes <- as.character(svm_pdata_classes)
svm_pdata_classes[svm_pdata_classes == "X1"] <- "1"
svm_pdata_classes[svm_pdata_classes == "X2"] <- "2"
svm_pdata_classes <- factor(svm_pdata_classes, levels = levels(testing_pdata$DX_GROUP))

roc_pdata <- roc(
  testing_pdata$DX_GROUP,
  svm_pdata_predictions[, 2],
  levels = rev(levels(testing_pdata$DX_GROUP))
)
auc_pdata <- auc(roc_pdata)
cat("\nAUC for testing_pdata SVM model:\n")
print(auc_pdata) #Area under the curve: 0.7659

conf_matrix_pdata <- confusionMatrix(svm_pdata_classes, testing_pdata$DX_GROUP)
cat("\nConfusion Matrix for testing_pdata SVM model:\n")
print(conf_matrix_pdata)

conf_matrix_pdata_df <- as.data.frame(as.table(conf_matrix_pdata$table))
conf_matrix_pdata_metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Value = c(conf_matrix_pdata$overall["Accuracy"],
            conf_matrix_pdata$byClass["Sensitivity"],
            conf_matrix_pdata$byClass["Specificity"])
)

conf_matrix_pdata_combined <- list(
  Confusion_Table = conf_matrix_pdata_df,
  Metrics = conf_matrix_pdata_metrics
)

conf_matrix_vdata_df <- as.data.frame(as.table(conf_matrix_vdata$table))
conf_matrix_vdata_metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Value = c(conf_matrix_vdata$overall["Accuracy"],
            conf_matrix_vdata$byClass["Sensitivity"],
            conf_matrix_vdata$byClass["Specificity"])
)

conf_matrix_vdata_combined <- list(
  Confusion_Table = conf_matrix_vdata_df,
  Metrics = conf_matrix_vdata_metrics
)

write.csv(conf_matrix_pdata_df, "conf_matrix_pdata_table.csv", row.names = FALSE)
write.csv(conf_matrix_pdata_metrics, "conf_matrix_pdata_metrics.csv", row.names = FALSE)
write.csv(conf_matrix_vdata_df, "conf_matrix_vdata_table.csv", row.names = FALSE)
write.csv(conf_matrix_vdata_metrics, "conf_matrix_vdata_metrics.csv", row.names = FALSE)

predictions_vdata <- as.data.frame(svm_vdata_classes)
predictions_pdata <- as.data.frame(svm_pdata_classes)
svm_oos_output_vdata <- cbind(testing_vdata, predictions_vdata)
svm_oos_output_pdata <- cbind(testing_pdata, predictions_pdata)

write.csv(svm_oos_output_vdata, "svm_oos_output_vdata.csv")
write.csv(svm_oos_output_pdata, "svm_oos_output_pdata.csv")

###IN-SAMPLE / OUT-OF-SAMPLE ROC CURVES
par(mfrow = c(2, 2))

roc_cv_vdata <- roc(
  response = svm_vdata$pred$obs,
  predictor = svm_vdata$pred$X1,
  levels = rev(levels(svm_vdata$pred$obs))
)
plot(roc_cv_vdata, col = "blue", main = "In-Sample ROC - Varimax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_cv_vdata), 4)), col = "blue", lwd = 2)

roc_cv_pdata <- roc(
  response = svm_pdata$pred$obs,
  predictor = svm_pdata$pred$X1,
  levels = rev(levels(svm_pdata$pred$obs))
)
plot(roc_cv_pdata, col = "darkgreen", main = "In-Sample ROC - Promax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_cv_pdata), 4)), col = "darkgreen", lwd = 2)

roc_oos_vdata <- roc(
  response = testing_vdata$DX_GROUP,
  predictor = svm_vdata_predictions[, 2],
  levels = rev(levels(testing_vdata$DX_GROUP))
)
plot(roc_oos_vdata, col = "red", main = "Out-of-Sample ROC - Varimax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_oos_vdata), 4)), col = "red", lwd = 2)

roc_oos_pdata <- roc(
  response = testing_pdata$DX_GROUP,
  predictor = svm_pdata_predictions[, 2],
  levels = rev(levels(testing_pdata$DX_GROUP))
)
plot(roc_oos_pdata, col = "purple", main = "Out-of-Sample ROC - Promax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_oos_pdata), 4)), col = "purple", lwd = 2)

})

## ================================================================================
## STAGE 21 — SVM MANCOVA (Actual DX_GROUP): Communication  (source: SVM_MANCOVA_comm_Actual.R)
## ================================================================================
run_stage(21, "SVM MANCOVA (Actual) - Communication", {

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "actual_all_output_detailed_vdf[comm].csv")
write.csv(varimax_summary_df, "actual_all_output_summary_vdf[comm].csv")
write.csv(promax_detailed_summary_df, "actual_all_output_detailed_pdf[comm].csv")
write.csv(promax_summary_df, "actual_all_output_summary_pdf[comm].csv")

})

## ================================================================================
## STAGE 22 — SVM MANCOVA (Model Predicted): Communication  (source: SVM_MANCOVA_comm_Model.R)
## FLAG: original renames the SVM prediction column to `pred` but the formula
## still references `DX_GROUP`, so this stage currently computes the exact same
## MANCOVA as Stage 21. Preserved as originally written; swap DX_GROUP -> pred
## in the formula below if the intent was to model against predicted labels.
## ================================================================================
run_stage(22, "SVM MANCOVA (Model) - Communication", {

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

colnames(vdata)[colnames(vdata) == "svm_vdata_classes_all"] <- "pred"
colnames(pdata)[colnames(pdata) == "svm_pdata_classes_all"] <- "pred"

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_SOCIAL_INT +
        SEX*DX_GROUP + SEX*AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "model_all_output_detailed_vdf[comm].csv")
write.csv(varimax_summary_df, "model_all_output_summary_vdf[comm].csv")
write.csv(promax_detailed_summary_df, "model_all_output_detailed_pdf[comm].csv")
write.csv(promax_summary_df, "model_all_output_summary_pdf[comm].csv")

})

## ================================================================================
## STAGE 23 — SVM MANCOVA (Actual DX_GROUP): Social  (source: SVM_MANCOVA_social_Actual.R)
## Note: original reads the out-of-sample-only files (svm_oos_output_vdata.csv /
## svm_oos_output_pdata.csv), not the "_all" versions used elsewhere — preserved.
## ================================================================================
run_stage(23, "SVM MANCOVA (Actual) - Social", {

vdata <- fread("svm_oos_output_vdata.csv")
pdata <- fread("svm_oos_output_pdata.csv")

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "actual_oos_output_detailed_vdf[social].csv")
write.csv(varimax_summary_df, "actual_oos_output_summary_vdf[social].csv")
write.csv(promax_detailed_summary_df, "actual_oos_output_detailed_pdf[social].csv")
write.csv(promax_summary_df, "actual_oos_output_summary_pdf[social].csv")

})

## ================================================================================
## STAGE 24 — SVM MANCOVA (Model Predicted): Social  (source: SVM_MANCOVA_social_Model.R)
## FLAG: same issue as Stage 22 — renames to `pred` but formula still uses
## `DX_GROUP`, so this duplicates Stage 21/23's logic. Preserved as written.
## ================================================================================
run_stage(24, "SVM MANCOVA (Model) - Social", {

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

colnames(vdata)[colnames(vdata) == "svm_vdata_classes_all"] <- "pred"
colnames(pdata)[colnames(pdata) == "svm_pdata_classes_all"] <- "pred"

dependent_vars <- paste0("PC", 1:20)

mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "),
        ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_SOCIAL_INT +
        SEX*DX_GROUP + SEX*AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
)

varimax_mancova_result <- manova(mancova_formula, data = vdata)
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")
varimax_summary_aov <- summary.aov(varimax_mancova_result)

promax_mancova_result <- manova(mancova_formula, data = pdata)
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL

varimax_anova_list <- lapply(varimax_summary_aov, function(x) as.data.frame(x))
promax_anova_list <- lapply(promax_summary_aov, function(x) as.data.frame(x))

varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var
  return(df)
}))

write.csv(varimax_detailed_summary_df, "model_all_output_detailed_vdf[social].csv")
write.csv(varimax_summary_df, "model_all_output_summary_vdf[social].csv")
write.csv(promax_detailed_summary_df, "model_all_output_detailed_pdf[social].csv")
write.csv(promax_summary_df, "model_all_output_summary_pdf[social].csv")

})

## ================================================================================
## STAGE 25 — SVM Permutation MANCOVA (Model Predicted): Communication  (source: SVM_PermMANCOVA_comm_Model.R)
## ================================================================================
run_stage(25, "SVM Permutation MANCOVA (Model) - Communication", {

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

colnames(vdata)[colnames(vdata) == "svm_vdata_classes_all"] <- "pred"
colnames(pdata)[colnames(pdata) == "svm_pdata_classes_all"] <- "pred"

pdata$pred <- as.factor(pdata$pred)
vdata$pred <- as.factor(vdata$pred)

if (length(levels(pdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in pdata.")
if (length(levels(vdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in vdata.")

dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, ..dv_names]
dependent_vars_v <- vdata[, ..dv_names]

set.seed(123)

perm_mancova_p <- anova.cca(capscale(as.matrix(dependent_vars_p) ~ pred * AVG_COMM_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_COMM_INT +
                                       SEX * pred + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                                     data = pdata), by = "margin", permutations = 10000)
print(perm_mancova_p)
write.csv(perm_mancova_p, "PermMANCOVA_svmdata_promax_comm.csv")

perm_mancova_v <- anova.cca(capscale(as.matrix(dependent_vars_v) ~ pred * AVG_COMM_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_COMM_INT +
                                       SEX * pred + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                                     data = vdata), by = "margin", permutations = 10000)
print(perm_mancova_v)
write.csv(perm_mancova_v, "PermMANCOVA_svmdata_varimax_comm.csv")

bonferroni_interaction_tests <- function(data, predictor, interaction_term) {
  results <- list()
  for (dv in dv_names) {
    formula <- as.formula(paste(dv, "~", predictor, "*", interaction_term))
    bonferroni_result <- tryCatch({
      pairwise.t.test(data[[dv]], interaction(data[[predictor]], data[[interaction_term]]), p.adjust.method = "bonferroni")
    }, error = function(e) {
      message(paste("Bonferroni test failed for", dv, "due to insufficient group levels"))
      return(NULL)
    })
    if (!is.null(bonferroni_result)) {
      results[[paste(dv, "Bonferroni Interaction")]] <- bonferroni_result$p.value
    }
  }
  return(results)
}

posthoc_p_bonferroni_interaction <- bonferroni_interaction_tests(pdata, "pred", "AVG_COMM_INT")
posthoc_v_bonferroni_interaction <- bonferroni_interaction_tests(vdata, "pred", "AVG_COMM_INT")

print(posthoc_p_bonferroni_interaction)
print(posthoc_v_bonferroni_interaction)

})

## ================================================================================
## STAGE 26 — SVM Permutation MANCOVA (Model Predicted): Social  (source: SVM_PermMANCOVA_social_Model.R)
## ================================================================================
run_stage(26, "SVM Permutation MANCOVA (Model) - Social", {

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

colnames(vdata)[colnames(vdata) == "svm_vdata_classes_all"] <- "pred"
colnames(pdata)[colnames(pdata) == "svm_pdata_classes_all"] <- "pred"

pdata$pred <- as.factor(pdata$pred)
vdata$pred <- as.factor(vdata$pred)

if (length(levels(pdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in pdata.")
if (length(levels(vdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in vdata.")

dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, ..dv_names]
dependent_vars_v <- vdata[, ..dv_names]

set.seed(123)

perm_mancova_p <- anova.cca(capscale(as.matrix(dependent_vars_p) ~ pred * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * pred + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = pdata), by = "margin", permutations = 10000)
print(perm_mancova_p)
write.csv(perm_mancova_p, "PermMANCOVA_svmdata_promax_social.csv")

perm_mancova_v <- anova.cca(capscale(as.matrix(dependent_vars_v) ~ pred * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * pred + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = vdata), by = "margin", permutations = 10000)
print(perm_mancova_v)
write.csv(perm_mancova_v, "PermMANCOVA_svmdata_varimax_social.csv")

bonferroni_interaction_tests <- function(data, predictor, interaction_term) {
  results <- list()
  for (dv in dv_names) {
    formula <- as.formula(paste(dv, "~", predictor, "*", interaction_term))
    bonferroni_result <- tryCatch({
      pairwise.t.test(data[[dv]], interaction(data[[predictor]], data[[interaction_term]]), p.adjust.method = "bonferroni")
    }, error = function(e) {
      message(paste("Bonferroni test failed for", dv, "due to insufficient group levels"))
      return(NULL)
    })
    if (!is.null(bonferroni_result)) {
      results[[paste(dv, "Bonferroni Interaction")]] <- bonferroni_result$p.value
    }
  }
  return(results)
}

posthoc_p_bonferroni_interaction <- bonferroni_interaction_tests(pdata, "pred", "AVG_SOCIAL_INT")
posthoc_v_bonferroni_interaction <- bonferroni_interaction_tests(vdata, "pred", "AVG_SOCIAL_INT")

print(posthoc_p_bonferroni_interaction)
print(posthoc_v_bonferroni_interaction)

})

## ================================================================================
## STAGE 27 — SVM Post-Hoc Tests: Tukey HSD & Games-Howell  (source: SVM_PermMANCOVA_Post Hoc.R)
## ================================================================================
run_stage(27, "SVM Post-Hoc Tests (Tukey HSD / Games-Howell)", {

pdata <- fread("svm_oos_output_pdata_all.csv")
vdata <- fread("svm_oos_output_vdata_all.csv")

pc_columns <- paste0("PC", 1:20)

pdata$pred <- as.factor(pdata$svm_pdata_classes_all)
vdata$pred <- as.factor(vdata$svm_vdata_classes_all)

#Tukey's HSD for AVG_SOCIAL_INT
pdata$AVG_SOCIAL_INT <- as.factor(pdata$AVG_SOCIAL_INT)
vdata$AVG_SOCIAL_INT <- as.factor(vdata$AVG_SOCIAL_INT)

run_tukey_hsd <- function(data, dataset_name, factor_var) {
  results <- list()

  for (dv in pc_columns) {
    if (!dv %in% colnames(data)) {
      message(paste("Skipping", dv, "because it is not found in the dataset."))
      next
    }

    formula <- as.formula(paste(dv, "~ pred *", factor_var))

    model <- tryCatch({
      aov(formula, data = data)
    }, error = function(e) {
      message(paste("ANOVA failed for", dv, ":", e$message))
      return(NULL)
    })

    if (!is.null(model)) {
      message(paste("Running Tukey's HSD for", dv))
      results[[dv]] <- tryCatch({
        TukeyHSD(model)
      }, error = function(e) {
        message(paste("Tukey's HSD failed for", dv, ":", e$message))
        return(NULL)
      })
    }
  }

  tidy_results <- list()
  for (dv in names(results)) {
    if (!is.null(results[[dv]])) {
      tidy_data <- as.data.frame(results[[dv]][[1]]) %>%
        tibble::rownames_to_column(var = "Comparison") %>%
        mutate(Dependent_Variable = dv)
      tidy_results[[dv]] <- tidy_data
    }
  }

  final_results <- bind_rows(tidy_results)
  final_results$Dataset <- dataset_name
  return(final_results)
}

tukey_results_pdata <- run_tukey_hsd(pdata, "pdata", "AVG_SOCIAL_INT")
tukey_results_vdata <- run_tukey_hsd(vdata, "vdata", "AVG_SOCIAL_INT")

combined_results <- bind_rows(tukey_results_pdata, tukey_results_vdata)
print(combined_results)
write.csv(combined_results, "Tukey_HSD_Results_Social.csv", row.names = FALSE)

#Tukey's HSD for AVG_COMM_INT
pdata$AVG_COMM_INT <- as.factor(pdata$AVG_COMM_INT)
vdata$AVG_COMM_INT <- as.factor(vdata$AVG_COMM_INT)

tukey_results_pdata <- run_tukey_hsd(pdata, "pdata", "AVG_COMM_INT")
tukey_results_vdata <- run_tukey_hsd(vdata, "vdata", "AVG_COMM_INT")

combined_results <- bind_rows(tukey_results_pdata, tukey_results_vdata)
print(combined_results)
write.csv(combined_results, "Tukey_HSD_Results_Comm.csv", row.names = FALSE)

run_games_howell <- function(data, factor_var, dataset_name) {
  results <- list()

  for (dv in pc_columns) {
    required_cols <- c(dv, "pred", factor_var)
    if (!all(required_cols %in% names(data))) {
      message(paste("Skipping", dv, "- missing required columns."))
      next
    }

    subset_data <- data[, ..required_cols]
    subset_data <- subset_data[complete.cases(subset_data), ]

    subset_data$group <- interaction(subset_data$pred, subset_data[[factor_var]], drop = TRUE)

    group_counts <- table(subset_data$group)
    valid_groups <- names(group_counts[group_counts >= 2])
    subset_data <- subset_data[subset_data$group %in% valid_groups, ]

    if (length(unique(subset_data$group)) >= 2) {
      gh_result <- tryCatch({
        oneway.test(subset_data[[dv]] ~ subset_data$group, var.equal = FALSE)
      }, error = function(e) {
        message(paste("Games-Howell test failed for", dv, ":", e$message))
        return(NULL)
      })

      if (!is.null(gh_result)) {
        results[[dv]] <- data.frame(
          Dependent_Variable = dv,
          F = gh_result$statistic,
          df = gh_result$parameter[2],
          p_value = gh_result$p.value
        )
      }
    } else {
      message(paste("Skipping", dv, "- not enough valid groups"))
    }
  }

  final_results <- dplyr::bind_rows(results)
  final_results$Dataset <- dataset_name
  return(final_results)
}

gh_pdata_social <- run_games_howell(pdata, "AVG_SOCIAL_INT", "pdata")
gh_vdata_social <- run_games_howell(vdata, "AVG_SOCIAL_INT", "vdata")

gh_pdata_comm <- run_games_howell(pdata, "AVG_COMM_INT", "pdata")
gh_vdata_comm <- run_games_howell(vdata, "AVG_COMM_INT", "vdata")

final_social <- bind_rows(gh_pdata_social, gh_vdata_social)
final_comm <- bind_rows(gh_pdata_comm, gh_vdata_comm)

write.csv(final_social, "Games_Howell_Results_Social.csv", row.names = FALSE)
write.csv(final_comm, "Games_Howell_Results_Comm.csv", row.names = FALSE)

})

message("\n===== PIPELINE COMPLETE =====")
