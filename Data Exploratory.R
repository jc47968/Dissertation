library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)


setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\05.1 Data Exporatory]")

# Load datasets
promax_data <- fread("promax_alldata.csv")
varimax_data <- fread("varimax_alldata.csv")

# Select PC1 to PC20 columns
pc_columns <- paste0("PC", 1:20)

# Convert DX_GROUP to factor
promax_data$DX_GROUP <- as.factor(promax_data$DX_GROUP)
varimax_data$DX_GROUP <- as.factor(varimax_data$DX_GROUP)

### QQ PLOT for PC1 to PC20 (Promax & Varimax) ###
# Function to generate QQ plots
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

# Generate and display QQ plots
qq_plot_promax <- generate_qq_plot(promax_data, "Promax")
qq_plot_varimax <- generate_qq_plot(varimax_data, "Varimax")

print(qq_plot_promax)
print(qq_plot_varimax)

### CDF PLOT for PC1 to PC20 by DX_GROUP ###
# Function to generate CDF plots
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

# Generate and display CDF plots
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

# Run KS test for both datasets
ks_results_promax <- perform_ks_test(promax_data, "Promax")
ks_results_varimax <- perform_ks_test(varimax_data, "Varimax")

# Combine results
ks_results <- rbind(ks_results_promax, ks_results_varimax)

# Print KS Test Results
print(ks_results)

# Save KS test results as CSV
write.csv(ks_results, "KS_Test_Results.csv", row.names = FALSE)
