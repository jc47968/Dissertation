library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(corrplot)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\EDA_Initial]")
data <- fread("merged_data.csv")
varimax_alldata <- fread("varimax_alldata.csv")
promax_alldata <- fread("promax_alldata.csv")
data <- data %>% select(SITE_ID:AVG_COMM)


# Summary of the dataset
summary(data)

# Create a table grouped by DX_GROUP for AGE_AT_SCAN and SEX variables
dx_group_table <- data %>%
  group_by(DX_GROUP) %>%
  summarize(
    Mean_Age = mean(AGE_AT_SCAN, na.rm = TRUE),
    SD_Age = sd(AGE_AT_SCAN, na.rm = TRUE),
    Male_Count = sum(SEX == 1, na.rm = TRUE),
    Female_Count = sum(SEX == 2, na.rm = TRUE)
  )
# Save the table to a CSV file
write.csv(dx_group_table, file = "dx_group_table.csv", row.names = FALSE)

# Distribution of Age
age_plot <- ggplot(data, aes(x = AGE_AT_SCAN)) +
  geom_histogram(binwidth = 2, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of Age", x = "Age", y = "Frequency")
print(age_plot)
ggsave("age_plot.png", plot = age_plot)

# Distribution of Diagnostic Group
dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP))) +
  geom_bar(fill = "green", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
  labs(title = "Distribution of Diagnostic Groups", x = "Diagnostic Group", y = "Count")
print(dx_plot)
ggsave("dx_plot.png", plot = dx_plot)

# Gender Distribution
gender_plot <- ggplot(data, aes(x = as.factor(SEX))) +
  geom_bar(fill = "purple", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")
print(gender_plot)
ggsave("gender_plot.png", plot = gender_plot)

# Pairplot to visualize relationships between key variables
pairplot_vars <- data %>% select(DX_GROUP, AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM)
png("pairplot.png")
pairs(pairplot_vars, main = "Pairplot of Selected Key Variables")
dev.off()

# Boxplot of Age by Diagnostic Group
age_dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP), y = AGE_AT_SCAN)) +
  geom_boxplot(fill = "orange", alpha = 0.7) +
  scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
  labs(title = "Age by Diagnostic Group", x = "Diagnostic Group", y = "Age")
print(age_dx_plot)
ggsave("age_dx_plot.png", plot = age_dx_plot)

# Relationship between AVG_SOCIAL and AGE_AT_SCAN
social_age_plot <- ggplot(data, aes(x = AGE_AT_SCAN, y = AVG_SOCIAL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "AVG_SOCIAL vs Age", x = "Age", y = "AVG_SOCIAL")
print(social_age_plot)
ggsave("social_age_plot.png", plot = social_age_plot)

# Visualization for AVG_SOCIAL and SEX
avg_social_sex_plot <- ggplot(data, aes(x = as.factor(SEX), y = AVG_SOCIAL)) +
  geom_boxplot(fill = "lightgreen", alpha = 0.7) +
  scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
  labs(title = "AVG_SOCIAL by Sex", x = "Sex", y = "AVG_SOCIAL")
print(avg_social_sex_plot)
ggsave("avg_social_sex_plot.png", plot = avg_social_sex_plot)

# Analysis of relationship between DX_GROUP and selected variables
selected_vars <- c('ADI_R_SOCIAL_TOTAL_A', 'ADI_R_VERBAL_TOTAL_BV', 'ADOS_COMM', 
                   'ADOS_SOCIAL', 'SRS_COGNITION', 'SRS_COMMUNICATION', 
                   'AVG_SOCIAL', 'AVG_COMM')

# Summary statistics by DX_GROUP
summary_stats <- data %>%
  group_by(DX_GROUP) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(summary_stats)

# Visualizations for each variable by DX_GROUP
for (var in selected_vars) {
  dx_plot <- ggplot(data, aes_string(x = "as.factor(DX_GROUP)", y = var)) +
    geom_boxplot(fill = "cyan", alpha = 0.7) +
    scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
    labs(title = paste(var, "by Diagnostic Group"), x = "Diagnostic Group", y = var)
  ggsave(paste0("dx_plot_", var, ".png"), plot = dx_plot)
}

# Pairwise comparison for DX_GROUP and selected variables using ANOVA
anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ DX_GROUP")), data = data)
  anova_results[[var]] <- summary(model)
}
print(anova_results)

dx_group_pvalues <- sapply(anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(dx_group_pvalues), file = "dx_group_pvalues.csv", row.names = TRUE)

# Individual variable trend visualizations
for (var in selected_vars) {
  trend_plot <- ggplot(data, aes_string(x = "AGE_AT_SCAN", y = var, color = "as.factor(DX_GROUP)")) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste(var, "Trends by Age and Diagnostic Group"), x = "Age", y = var, color = "Diagnostic Group")
  print(trend_plot)
}

# Analysis of relationship between SEX and selected variables
# Summary statistics by SEX
sex_summary_stats <- data %>%
  group_by(SEX) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(sex_summary_stats)

# Visualizations for each variable by SEX
for (var in selected_vars) {
  sex_plot <- ggplot(data, aes_string(x = "as.factor(SEX)", y = var)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
    labs(title = paste(var, "by Sex"), x = "Sex", y = var)
  ggsave(paste0("sex_plot_", var, ".png"), plot = sex_plot)
}

# Pairwise comparison for SEX and selected variables using ANOVA
sex_anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ SEX")), data = data)
  sex_anova_results[[var]] <- summary(model)
}
print(sex_anova_results)
# Extract and save p-values for sex_anova_results results
sex_anova_pvalues <- sapply(sex_anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(sex_anova_pvalues), file = "sex_anova_pvalues.csv", row.names = TRUE)

# Analysis of relationship between AGE_AT_SCAN and selected variables
# Summary statistics by AGE_AT_SCAN
data <- data %>% mutate(AGE_INTERVAL = cut(AGE_AT_SCAN, 
                                           breaks = c(-Inf, 10, 20, 30, 40, Inf), 
                                           labels = c("<10", "10-20", "20-30", "30-40", ">40")))

AGE_INTERVAL_summary_stats <- data %>%
  group_by(AGE_INTERVAL) %>%
  summarize(across(all_of(selected_vars), list(mean = mean, sd = sd), na.rm = TRUE))
print(AGE_INTERVAL_summary_stats)

# Visualizations for each variable by AGE_INTERVAL
for (var in selected_vars) {
  AGE_INTERVAL_plot <- ggplot(data, aes(x = AGE_INTERVAL, y = .data[[var]])) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    labs(title = paste(var, "by AGE_INTERVAL"), x = "AGE_INTERVAL", y = var)
  ggsave(paste0("AGE_INTERVAL_plot_", var, ".png"), plot = AGE_INTERVAL_plot)
}

# Pairwise comparison for AGE_INTERVAL and selected variables using ANOVA
AGE_INTERVAL_anova_results <- list()
for (var in selected_vars) {
  model <- aov(as.formula(paste(var, "~ AGE_INTERVAL")), data = data)
  AGE_INTERVAL_anova_results[[var]] <- summary(model)
}
print(AGE_INTERVAL_anova_results)
# Extract and save p-values for AGE_INTERVAL ANOVA results
age_interval_pvalues <- sapply(AGE_INTERVAL_anova_results, function(x) x[[1]]["Pr(>F)"][1])
write.csv(as.data.frame(age_interval_pvalues), file = "age_interval_anova_pvalues.csv", row.names = TRUE)

#  AGE_INTERVAL for Social and Communication Variables
age_int_grouped_table <- data %>%
  group_by(AGE_INTERVAL) %>%
  summarize(across(ADI_R_SOCIAL_TOTAL_A:AVG_COMM, mean, na.rm = TRUE))

# Save the table to a CSV file
write.csv(age_int_grouped_table, file = "grouped_table_by_age_interval.csv", row.names = FALSE)

# Select PC1 to PC20 from varimax_alldata and promax_alldata
varimax_data <- varimax_alldata %>% select(where(is.numeric)) %>% select(PC1:PC20)
promax_data <- promax_alldata %>% select(where(is.numeric)) %>% select(PC1:PC20)

# Compute correlation matrices
varimax_corr <- cor(varimax_data, use = "complete.obs")
promax_corr <- cor(promax_data, use = "complete.obs")

# Plot correlation matrix for varimax_alldata
png("varimax_correlation_matrix.png", width = 800, height = 800)
corrplot(varimax_corr, method = "color", title = "Correlation Matrix: Varimax", tl.cex = 0.8, cl.cex = 0.8, addgrid.col = NA)
dev.off()

# Plot correlation matrix for promax_alldata
png("promax_correlation_matrix.png", width = 800, height = 800)
corrplot(promax_corr, method = "color", title = "Correlation Matrix: Promax", tl.cex = 0.8, cl.cex = 0.8, addgrid.col = NA)
dev.off()

##PC Average for DX_GROUP
# Load the datasets
promax_alldata <- read.csv("promax_alldata.csv")
varimax_alldata <- read.csv("varimax_alldata.csv")

# Define a function to compute summary statistics for PC columns grouped by DX_GROUP
summarize_pcs_by_group <- function(data, dataset_name) {
  data %>%
    # Select PC1 to PC20 and DX_GROUP
    select(DX_GROUP, PC1:PC20) %>%
    # Group by DX_GROUP
    group_by(DX_GROUP) %>%
    # Compute mean for each PC within each group
    summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
    # Pivot longer for PC rows
    pivot_longer(-DX_GROUP, names_to = "PC", values_to = "Mean") %>%
    # Pivot wider for DX_GROUP columns
    pivot_wider(names_from = DX_GROUP, values_from = Mean) %>%
    # Rename columns for clarity
    rename_with(~ paste0(dataset_name, "_", .), -PC)
}

# Summarize PCs for promax_alldata
promax_summary <- summarize_pcs_by_group(promax_alldata, "Promax")

# Summarize PCs for varimax_alldata
varimax_summary <- summarize_pcs_by_group(varimax_alldata, "Varimax")

# Combine the results from both datasets
combined_summary <- full_join(promax_summary, varimax_summary, by = "PC")

write.csv(combined_summary, file = "pc_summary_by_dx_group.csv", row.names = FALSE)


# Prepare data for visualization
visualization_data <- combined_summary %>%
  pivot_longer(cols = -PC, names_to = "Dataset_DX_GROUP", values_to = "Mean") %>%
  separate(Dataset_DX_GROUP, into = c("Dataset", "DX_GROUP"), sep = "_") %>%
  mutate(DX_GROUP = ifelse(DX_GROUP == "1", "Autism", "Control")) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:20)))

# Create a grouped bar plot
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

# Prepare data for CDF visualization
prepare_cdf_data <- function(data, dataset_name) {
  data %>%
    select(starts_with("PC")) %>%
    pivot_longer(cols = everything(), names_to = "PC", values_to = "Value") %>%
    filter(!is.na(Value) & PC %in% paste0("PC", 1:20)) %>% # Exclude NA values and ensure only PC1 to PC20
    group_by(PC) %>%
    mutate(CDF = ecdf(Value)(Value)) %>%
    ungroup() %>%
    mutate(Dataset = dataset_name)
}

# Prepare CDF data for promax_alldata and varimax_alldata
promax_cdf_data <- prepare_cdf_data(promax_alldata, "Promax")
varimax_cdf_data <- prepare_cdf_data(varimax_alldata, "Varimax")

# Combine both datasets
cdf_data <- bind_rows(promax_cdf_data, varimax_cdf_data) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:20)))

# Perform KS test for each PC
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

# Combine KS test results into a single data frame
ks_test_summary <- do.call(rbind, ks_test_results)

# Save KS test results to a CSV file
write.csv(ks_test_summary, file = "ks_test_results.csv", row.names = FALSE)

print(ks_test_summary)

# Create the CDF plot
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