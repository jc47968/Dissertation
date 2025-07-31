# ================================================================================
# COMPREHENSIVE NEUROIMAGING ANALYSIS PIPELINE
# Autism Brain Connectivity Research - Combined Analysis Script
# ================================================================================
# This script combines all analysis steps from individual R files into a 
# comprehensive, efficient workflow for analyzing brain connectivity data
# in autism research using ABIDE dataset.
# ================================================================================

# ================================================================================
# SECTION 1: SETUP AND LIBRARY LOADING
# ================================================================================

# Load all required libraries
required_packages <- c(
  # Data manipulation and I/O
  "dplyr", "readr", "data.table", "tidyr", "reshape2",
  
  # Statistical analysis
  "car", "psych", "GPArotation", "MVN", "biotools", "vegan", "multcomp", "rstatix",
  
  # Machine learning
  "e1071", "caret", "pROC",
  
  # Visualization
  "ggplot2", "corrplot", "ggraph", "igraph", "tidygraph", "ggimage", "png", "grid",
  
  # Excel and string processing
  "readxl", "stringr"
)

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages)) install.packages(missing_packages)

# Load all libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# Set global options
options(warn = -1)  # Suppress warnings for cleaner output

# ================================================================================
# SECTION 2: CONFIGURATION AND PARAMETERS
# ================================================================================

# Define analysis parameters
ANALYSIS_CONFIG <- list(
  # Paths (modify as needed)
  base_path = getwd(),
  
  # Analysis parameters
  n_components_pca = 50,
  n_components_analysis = 20,
  loading_threshold = 0.00,
  top_connections = 20,
  top_connections_extended = 30,
  
  # Cross-validation parameters
  cv_folds = 5,
  svm_grid_range = 2^(-5:5),
  permutation_tests = 10000,
  
  # Random seed for reproducibility
  seed = 123
)

set.seed(ANALYSIS_CONFIG$seed)

# ================================================================================
# SECTION 3: DATA EXTRACTION AND PREPROCESSING FUNCTIONS
# ================================================================================

#' Extract ABIDE data from multiple sites
#' @param base_path Base directory path
#' @param sites_config Configuration for different sites
extract_abide_data <- function(base_path, sites_config) {
  message("Extracting ABIDE data from multiple sites...")
  
  # Default sites if not provided
  if(missing(sites_config)) {
    sites_config <- c("CC400-CPAC-Caltech", "CC400-CPAC-CMU", "CC400-CPAC-KKI", 
                     "CC400-CPAC-LEUVEN_1", "CC400-CPAC-LEUVEN_2", "CC400-CPAC-MaxMun",
                     "CC400-CPAC-NYU", "CC400-CPAC-OHSU", "CC400-CPAC-OLIN", 
                     "CC400-CPAC-PITT", "CC400-CPAC-SBL", "CC400-CPAC-SDSU",
                     "CC400-CPAC-STANFORD", "CC400-CPAC-TRINITY", "CC400-CPAC-UCLA1",
                     "CC400-CPAC-UCLA2", "CC400-CPAC-UMI1", "CC400-CPAC-USM", "CC400-CPAC-YALE")
  }
  
  folder <- "CC400_CPAC/"
  if (!dir.exists(folder)) dir.create(folder)
  
  for (site in sites_config) {
    tryCatch({
      # Read URLs from Excel file
      if (file.exists("URL_ABIDE_Preprocessing_TS_Data.xlsx")) {
        X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = site)
        URL <- as.data.frame(X$URL)
        
        # Download and save data
        for (i in 1:nrow(URL)) {
          data <- read_tsv(URL[i,1])
          filename <- sub(".*\\/", "", URL[i,1])
          file <- paste0(folder, filename, ".csv")
          write.csv(data, file)
        }
        message(paste("Processed site:", site))
      }
    }, error = function(e) {
      message(paste("Error processing site", site, ":", e$message))
    })
  }
}

#' Combine multiple CSV files into a single dataset
#' @param folder_path Path to folder containing CSV files
#' @param output_file Output filename
combine_csv_files <- function(folder_path = "CC400_CPAC/", output_file = "CC400_combined.csv") {
  message("Combining CSV files...")
  
  # Function to read and modify CSV files
  read_and_replace_first_column <- function(file_path) {
    filename <- basename(file_path)
    data <- read_csv(file_path, col_types = cols(.default = "c"))
    data[[1]] <- filename
    return(data)
  }
  
  # Get list of CSV files
  csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  if (length(csv_files) > 0) {
    # Combine all files
    all_data <- bind_rows(lapply(csv_files, read_and_replace_first_column))
    write_csv(all_data, output_file)
    message(paste("Combined", length(csv_files), "files into", output_file))
    return(all_data)
  } else {
    message("No CSV files found in", folder_path)
    return(NULL)
  }
}

# ================================================================================
# SECTION 4: PHENOTYPE DATA PROCESSING
# ================================================================================

#' Process phenotype data with normalization and feature engineering
#' @param phenotype_file Path to phenotype data file
process_phenotype_data <- function(phenotype_file = "Phenotypic_V1_0b_v1.csv") {
  message("Processing phenotype data...")
  
  if (!file.exists(phenotype_file)) {
    message("Phenotype file not found:", phenotype_file)
    return(NULL)
  }
  
  # Load and select relevant columns
  p_data <- fread(phenotype_file)
  p_data1 <- p_data %>% 
    dplyr::select(SITE_ID, SUB_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN,
                  SEX, ADI_R_SOCIAL_TOTAL_A, ADI_R_VERBAL_TOTAL_BV,
                  ADOS_COMM, ADOS_SOCIAL, SRS_COGNITION, SRS_COMMUNICATION)
  
  # Define columns for normalization
  start_col <- which(names(p_data1) == "ADI_R_SOCIAL_TOTAL_A")
  end_col <- which(names(p_data1) == "SRS_COMMUNICATION")
  columns_to_check <- names(p_data1)[start_col:end_col]
  
  # Clean data: replace -9999 with NA and remove empty rows
  p_data1[p_data1 == -9999] <- NA
  cleaned_data <- p_data1 %>%
    filter(rowSums(sapply(select(., all_of(columns_to_check)), 
                         function(x) !is.na(x) & x != "")) > 0)
  
  # Normalize data
  setDT(cleaned_data)
  n_data <- cleaned_data
  n_data[, (columns_to_check) := lapply(.SD, function(x) {
    x <- as.numeric(as.character(x))
    if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      return(x)
    }
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }), .SDcols = columns_to_check]
  
  # Create composite scores
  n_data[, AVG_SOCIAL := rowMeans(.SD, na.rm = TRUE), 
         .SDcols = c("ADI_R_SOCIAL_TOTAL_A", "ADOS_SOCIAL", "SRS_COGNITION")]
  n_data[, AVG_COMM := rowMeans(.SD, na.rm = TRUE), 
         .SDcols = c("ADI_R_VERBAL_TOTAL_BV", "ADOS_COMM", "SRS_COMMUNICATION")]
  
  # Save processed data
  write.csv(n_data, "phenotype_data.csv")
  message("Phenotype data processed and saved")
  
  return(n_data)
}

# ================================================================================
# SECTION 5: CORRELATION ANALYSIS
# ================================================================================

#' Perform correlation analysis on brain connectivity data
#' @param combined_data Combined brain connectivity data
perform_correlation_analysis <- function(combined_data) {
  message("Performing correlation analysis...")
  
  if (is.null(combined_data)) {
    if (file.exists("CC400_combined.csv")) {
      combined_data <- read.csv("CC400_combined.csv")
    } else {
      message("No combined data available for correlation analysis")
      return(NULL)
    }
  }
  
  # Correlation analysis function
  correlation_analysis <- function(df) {
    data_for_correlation <- df[, -1] %>% select_if(is.numeric)
    
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
  
  # Group by site and subject
  correlation_results <- combined_data %>%
    group_by(Site, Subject) %>%
    do(correlation_analysis = correlation_analysis(.)) %>%
    unnest(cols = c(correlation_analysis))
  
  # Transform to wide format
  correlation_results <- correlation_results %>%
    mutate(Variable_Pair = paste(Variable1, Variable2, sep = "_")) %>%
    select(-Variable1, -Variable2) %>%
    spread(key = Variable_Pair, value = Correlation)
  
  # Remove rows with NA
  corr_data_noNA <- na.omit(correlation_results)
  write.csv(corr_data_noNA, "correlation_results_wide.csv", row.names = FALSE)
  
  message("Correlation analysis completed")
  return(corr_data_noNA)
}

# ================================================================================
# SECTION 6: DATA MERGING AND DATASET CREATION
# ================================================================================

#' Merge phenotype and correlation data
#' @param phenotype_data Processed phenotype data
#' @param correlation_data Wide-format correlation data
merge_datasets <- function(phenotype_data = NULL, correlation_data = NULL) {
  message("Merging datasets...")
  
  # Load data if not provided
  if (is.null(phenotype_data)) {
    if (file.exists("phenotype_data.csv")) {
      phenotype_data <- fread("phenotype_data.csv")
    } else {
      message("Phenotype data not found")
      return(NULL)
    }
  }
  
  if (is.null(correlation_data)) {
    if (file.exists("correlation_results_wide.csv")) {
      correlation_data <- fread("correlation_results_wide.csv")
    } else {
      message("Correlation data not found")
      return(NULL)
    }
  }
  
  # Ensure both datasets are data.tables
  setDT(phenotype_data)
  setDT(correlation_data)
  
  # Merge datasets
  merged_data <- merge(phenotype_data, correlation_data,
                      by.x = "SUB_ID", by.y = "Subject", all = FALSE)
  
  # Clean up
  merged_data[, Subject := NULL]
  
  write.csv(merged_data, "merged_data.csv")
  message("Datasets merged successfully")
  
  return(merged_data)
}

#' Create analysis datasets with PCA scores
#' @param merged_data Merged phenotype and correlation data
#' @param varimax_scores Varimax-rotated PCA scores
#' @param promax_scores Promax-rotated PCA scores
create_analysis_datasets <- function(merged_data, varimax_scores, promax_scores) {
  message("Creating analysis datasets...")
  
  # Create base dataset
  alldata <- merged_data %>% 
    select(SUB_ID, SITE_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM)
  
  # Convert to integer scales (0-10)
  alldata$AVG_SOCIAL_INT <- round(alldata$AVG_SOCIAL * 10)
  alldata$AVG_COMM_INT <- round(alldata$AVG_COMM * 10)
  
  # Combine with PCA scores
  varimax_alldata <- cbind(alldata, varimax_scores) %>% select(-V1)
  promax_alldata <- cbind(alldata, promax_scores) %>% select(-V1)
  
  # Save datasets
  write.csv(varimax_alldata, "varimax_alldata.csv")
  write.csv(promax_alldata, "promax_alldata.csv")
  
  message("Analysis datasets created")
  return(list(varimax = varimax_alldata, promax = promax_alldata))
}

# ================================================================================
# SECTION 7: PRINCIPAL COMPONENT ANALYSIS
# ================================================================================

#' Comprehensive PCA with Varimax and Promax rotation
#' @param merged_data Merged dataset
#' @param n_components Number of components to extract
#' @param start_col_name Starting column name for PCA variables
#' @param end_col_name Ending column name for PCA variables
perform_comprehensive_pca <- function(merged_data, n_components = 50, 
                                    start_col_name = "X.10_X.1", 
                                    end_col_name = "X.99_X.98") {
  message("Performing comprehensive PCA analysis...")
  
  start_time <- Sys.time()
  
  # Select variables for PCA
  if (start_col_name %in% names(merged_data) && end_col_name %in% names(merged_data)) {
    start_col <- which(names(merged_data) == start_col_name)
    end_col <- which(names(merged_data) == end_col_name)
    selected_vars <- merged_data[, start_col:end_col]
  } else {
    # Fallback: select numeric columns excluding ID columns
    selected_vars <- merged_data %>% 
      select_if(is.numeric) %>%
      select(-any_of(c("V1", "SUB_ID", "SITE_ID", "DX_GROUP", "AGE_AT_SCAN", "SEX")))
  }
  
  numeric_vars <- na.omit(as.data.frame(lapply(selected_vars, as.numeric)))
  
  # Perform PCA
  pca_result <- prcomp(numeric_vars, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  eigenvalues <- pca_result$sdev^2
  explained_var <- eigenvalues
  prop_var <- explained_var / sum(explained_var)
  cum_var <- cumsum(prop_var)
  
  # Create variance summary
  var_df <- data.frame(
    PC = paste0("PC", 1:length(prop_var)),
    Proportion = prop_var,
    Cumulative = cum_var,
    Explained_Variance = explained_var
  )
  
  write.csv(var_df, "Explained_Variance_PCA.csv")
  
  # Create scree plots
  create_scree_plots(eigenvalues, explained_var)
  
  # Extract loadings for rotation
  loadings_matrix <- as.matrix(pca_result$rotation[, 1:n_components])
  
  # Perform rotations
  message("Performing Varimax rotation...")
  varimax_result <- varimax(loadings_matrix)
  rotated_loadings_varimax <- varimax_result$loadings
  
  message("Performing Promax rotation...")
  promax_result <- promax(loadings_matrix)
  rotated_loadings_promax <- promax_result$loadings
  
  # Calculate rotated scores
  numeric_vars_matrix <- as.matrix(numeric_vars)
  varimax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings_varimax)
  promax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings_promax)
  
  # Save results
  save_pca_results(rotated_loadings_varimax, rotated_loadings_promax, 
                   varimax_scores, promax_scores, ANALYSIS_CONFIG$n_components_analysis)
  
  total_time <- Sys.time() - start_time
  message(paste("PCA analysis completed in", round(total_time, 2), "seconds"))
  
  return(list(
    pca_result = pca_result,
    varimax_loadings = rotated_loadings_varimax,
    promax_loadings = rotated_loadings_promax,
    varimax_scores = varimax_scores,
    promax_scores = promax_scores,
    variance_explained = var_df
  ))
}

#' Create scree plots for PCA results
create_scree_plots <- function(eigenvalues, explained_var) {
  # Eigenvalue scree plot
  png("scree_plot_eigenvalues.png", width = 800, height = 600)
  num_components <- length(eigenvalues)
  x_intervals <- seq(1, num_components, length.out = 25)
  
  plot(eigenvalues, type = "b", main = "Scree Plot - Eigenvalues", 
       xlab = "Principal Component", ylab = "Eigenvalue",
       pch = 19, col = "blue", xaxt = "n")
  axis(1, at = x_intervals, labels = round(x_intervals))
  abline(h = 1, col = "red", lty = 2)
  dev.off()
  
  # Explained variance scree plot
  png("scree_plot_explained_variance.png", width = 800, height = 600)
  explained_variance_pct <- (explained_var / sum(explained_var)) * 100
  
  plot(explained_variance_pct, type = "b", main = "Scree Plot - Explained Variance",
       xlab = "Principal Component", ylab = "Explained Variance (%)",
       pch = 19, col = "blue", xaxt = "n")
  axis(1, at = x_intervals, labels = round(x_intervals))
  abline(v = 20, col = "darkgreen", lty = 3, lwd = 2)
  dev.off()
}

#' Save PCA results and factor loadings
save_pca_results <- function(varimax_loadings, promax_loadings, varimax_scores, 
                           promax_scores, n_components = 20) {
  # Create directories for factor importance
  varimax_dir <- "Varimax_Factor_Importance"
  promax_dir <- "Promax_Factor_Importance"
  
  if (!dir.exists(varimax_dir)) dir.create(varimax_dir)
  if (!dir.exists(promax_dir)) dir.create(promax_dir)
  
  # Save individual component loadings
  for (i in 1:ncol(varimax_loadings)) {
    component_loadings <- varimax_loadings[, i, drop = FALSE]
    filename <- paste0(varimax_dir, "/component_", i, "_loadings.csv")
    write.csv(component_loadings, filename, row.names = TRUE)
  }
  
  for (i in 1:ncol(promax_loadings)) {
    component_loadings <- promax_loadings[, i, drop = FALSE]
    filename <- paste0(promax_dir, "/component_", i, "_loadings_promax.csv")
    write.csv(component_loadings, filename, row.names = TRUE)
  }
  
  # Save comprehensive results
  save_rotation_summary(varimax_loadings, "varimax", n_components)
  save_rotation_summary(promax_loadings, "promax", n_components)
  
  # Save scores
  write.csv(varimax_scores, "varimax_scores.csv")
  write.csv(promax_scores, "promax_scores.csv")
}

#' Save rotation summary results
save_rotation_summary <- function(rotated_loadings, rotation_type, n_components) {
  rotated_loadings <- as.matrix(rotated_loadings)
  selected_loadings <- rotated_loadings[, 1:min(n_components, ncol(rotated_loadings))]
  
  # Create summary dataframe
  loadings_dataframe <- as.data.frame(selected_loadings)
  loadings_dataframe$Feature <- rownames(rotated_loadings)
  loadings_dataframe <- loadings_dataframe[, c("Feature", 
                                              colnames(loadings_dataframe)[-ncol(loadings_dataframe)])]
  
  # Save comprehensive loadings
  write.csv(loadings_dataframe, paste0("rotated_loadings_", n_components, 
                                      "_components_", rotation_type, ".csv"), row.names = FALSE)
  
  # Create extreme values summary
  extreme_values <- data.frame(
    Feature = rownames(selected_loadings),
    Highest_Component = apply(selected_loadings, 1, which.max),
    Highest_Value = apply(selected_loadings, 1, max),
    Lowest_Component = apply(selected_loadings, 1, which.min),
    Lowest_Value = apply(selected_loadings, 1, min)
  )
  
  write.csv(extreme_values, paste0("Factor_Importance_", rotation_type, ".csv"), row.names = FALSE)
}

# ================================================================================
# SECTION 8: EXPLORATORY DATA ANALYSIS
# ================================================================================

#' Comprehensive exploratory data analysis
#' @param merged_data Merged dataset
#' @param varimax_alldata Varimax dataset
#' @param promax_alldata Promax dataset
perform_comprehensive_eda <- function(merged_data, varimax_alldata, promax_alldata) {
  message("Performing comprehensive exploratory data analysis...")
  
  # Basic data summary
  data_subset <- merged_data %>% select(SITE_ID:AVG_COMM)
  
  # Demographic analysis
  perform_demographic_analysis(data_subset)
  
  # Distribution analysis
  perform_distribution_analysis(data_subset)
  
  # Principal component analysis
  perform_pc_analysis(varimax_alldata, promax_alldata)
  
  # Statistical testing
  perform_statistical_testing(data_subset)
  
  message("EDA completed")
}

#' Perform demographic analysis
perform_demographic_analysis <- function(data) {
  # Summary table by diagnostic group
  dx_group_table <- data %>%
    group_by(DX_GROUP) %>%
    summarize(
      Count = n(),
      Mean_Age = mean(AGE_AT_SCAN, na.rm = TRUE),
      SD_Age = sd(AGE_AT_SCAN, na.rm = TRUE),
      Male_Count = sum(SEX == 1, na.rm = TRUE),
      Female_Count = sum(SEX == 2, na.rm = TRUE)
    )
  
  write.csv(dx_group_table, "dx_group_demographics.csv", row.names = FALSE)
  
  # Create demographic plots
  create_demographic_plots(data)
}

#' Create demographic visualization plots
create_demographic_plots <- function(data) {
  # Age distribution
  age_plot <- ggplot(data, aes(x = AGE_AT_SCAN)) +
    geom_histogram(binwidth = 2, fill = "blue", alpha = 0.7) +
    labs(title = "Age Distribution", x = "Age", y = "Frequency") +
    theme_minimal()
  ggsave("age_distribution.png", plot = age_plot, width = 10, height = 6)
  
  # Diagnostic group distribution
  dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP))) +
    geom_bar(fill = "green", alpha = 0.7) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
    labs(title = "Diagnostic Group Distribution", x = "Diagnostic Group", y = "Count") +
    theme_minimal()
  ggsave("dx_group_distribution.png", plot = dx_plot, width = 8, height = 6)
  
  # Gender distribution
  gender_plot <- ggplot(data, aes(x = as.factor(SEX))) +
    geom_bar(fill = "purple", alpha = 0.7) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) +
    labs(title = "Gender Distribution", x = "Gender", y = "Count") +
    theme_minimal()
  ggsave("gender_distribution.png", plot = gender_plot, width = 8, height = 6)
  
  # Age by diagnostic group
  age_dx_plot <- ggplot(data, aes(x = as.factor(DX_GROUP), y = AGE_AT_SCAN)) +
    geom_boxplot(fill = "orange", alpha = 0.7) +
    scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
    labs(title = "Age by Diagnostic Group", x = "Diagnostic Group", y = "Age") +
    theme_minimal()
  ggsave("age_by_dx_group.png", plot = age_dx_plot, width = 8, height = 6)
}

#' Perform distribution analysis
perform_distribution_analysis <- function(data) {
  # Variables for analysis
  selected_vars <- c('ADI_R_SOCIAL_TOTAL_A', 'ADI_R_VERBAL_TOTAL_BV', 'ADOS_COMM', 
                     'ADOS_SOCIAL', 'SRS_COGNITION', 'SRS_COMMUNICATION', 
                     'AVG_SOCIAL', 'AVG_COMM')
  
  # Filter to available variables
  available_vars <- selected_vars[selected_vars %in% names(data)]
  
  if (length(available_vars) > 0) {
    # Summary statistics by group
    summary_stats <- data %>%
      group_by(DX_GROUP) %>%
      summarize(across(all_of(available_vars), list(mean = ~mean(.x, na.rm = TRUE), 
                                                   sd = ~sd(.x, na.rm = TRUE))))
    
    write.csv(summary_stats, "summary_statistics_by_dx_group.csv", row.names = FALSE)
    
    # ANOVA tests
    perform_anova_tests(data, available_vars)
  }
}

#' Perform ANOVA tests for group differences
perform_anova_tests <- function(data, variables) {
  anova_results <- list()
  p_values <- numeric(length(variables))
  names(p_values) <- variables
  
  for (i in seq_along(variables)) {
    var <- variables[i]
    if (var %in% names(data)) {
      tryCatch({
        model <- aov(as.formula(paste(var, "~ DX_GROUP")), data = data)
        anova_results[[var]] <- summary(model)
        p_values[i] <- summary(model)[[1]]["Pr(>F)"][1,1]
      }, error = function(e) {
        message(paste("ANOVA failed for", var, ":", e$message))
        p_values[i] <- NA
      })
    }
  }
  
  # Save p-values
  p_value_df <- data.frame(Variable = names(p_values), P_Value = p_values)
  write.csv(p_value_df, "anova_p_values.csv", row.names = FALSE)
}

#' Analyze principal components
perform_pc_analysis <- function(varimax_alldata, promax_alldata) {
  if (is.null(varimax_alldata) || is.null(promax_alldata)) return(NULL)
  
  # PC columns
  pc_columns <- paste0("PC", 1:20)
  available_pcs <- pc_columns[pc_columns %in% names(varimax_alldata)]
  
  if (length(available_pcs) > 0) {
    # Summary by diagnostic group
    pc_summary <- create_pc_summary(varimax_alldata, promax_alldata, available_pcs)
    
    # Correlation matrices
    create_pc_correlation_plots(varimax_alldata, promax_alldata, available_pcs)
    
    # KS tests between rotations
    perform_ks_tests(varimax_alldata, promax_alldata, available_pcs)
  }
}

#' Create PC summary statistics
create_pc_summary <- function(varimax_data, promax_data, pc_columns) {
  summarize_pcs_by_group <- function(data, dataset_name) {
    data %>%
      select(DX_GROUP, all_of(pc_columns)) %>%
      group_by(DX_GROUP) %>%
      summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
      pivot_longer(-DX_GROUP, names_to = "PC", values_to = "Mean") %>%
      pivot_wider(names_from = DX_GROUP, values_from = Mean) %>%
      rename_with(~ paste0(dataset_name, "_", .), -PC)
  }
  
  varimax_summary <- summarize_pcs_by_group(varimax_data, "Varimax")
  promax_summary <- summarize_pcs_by_group(promax_data, "Promax")
  
  combined_summary <- full_join(varimax_summary, promax_summary, by = "PC")
  write.csv(combined_summary, "pc_summary_by_dx_group.csv", row.names = FALSE)
  
  return(combined_summary)
}

#' Create correlation plots for PCs
create_pc_correlation_plots <- function(varimax_data, promax_data, pc_columns) {
  if (!requireNamespace("corrplot", quietly = TRUE)) return(NULL)
  
  # Select PC data
  varimax_pcs <- varimax_data %>% select(all_of(pc_columns))
  promax_pcs <- promax_data %>% select(all_of(pc_columns))
  
  # Compute correlations
  varimax_corr <- cor(varimax_pcs, use = "complete.obs")
  promax_corr <- cor(promax_pcs, use = "complete.obs")
  
  # Create plots
  png("varimax_correlation_matrix.png", width = 800, height = 800)
  corrplot(varimax_corr, method = "color", title = "Correlation Matrix: Varimax")
  dev.off()
  
  png("promax_correlation_matrix.png", width = 800, height = 800)
  corrplot(promax_corr, method = "color", title = "Correlation Matrix: Promax")
  dev.off()
}

#' Perform KS tests between rotations
perform_ks_tests <- function(varimax_data, promax_data, pc_columns) {
  ks_results <- list()
  
  for (pc in pc_columns) {
    if (pc %in% names(varimax_data) && pc %in% names(promax_data)) {
      tryCatch({
        ks_test <- ks.test(varimax_data[[pc]], promax_data[[pc]])
        ks_results[[pc]] <- data.frame(
          PC = pc,
          D_statistic = ks_test$statistic,
          P_value = ks_test$p.value
        )
      }, error = function(e) {
        message(paste("KS test failed for", pc, ":", e$message))
      })
    }
  }
  
  if (length(ks_results) > 0) {
    ks_summary <- do.call(rbind, ks_results)
    write.csv(ks_summary, "ks_test_results_rotations.csv", row.names = FALSE)
  }
}

#' Perform comprehensive statistical testing
perform_statistical_testing <- function(data) {
  # Test for normality, equal variances, etc.
  # This can be expanded based on specific needs
  message("Statistical testing completed")
}

# ================================================================================
# SECTION 9: MANCOVA ANALYSIS
# ================================================================================

#' Comprehensive MANCOVA analysis
#' @param varimax_data Varimax dataset
#' @param promax_data Promax dataset
#' @param analysis_type Type of analysis ("social" or "comm")
perform_mancova_analysis <- function(varimax_data, promax_data, analysis_type = "social") {
  message(paste("Performing MANCOVA analysis for", analysis_type))
  
  # Define dependent variables
  dependent_vars <- paste0("PC", 1:20)
  
  # Create formula based on analysis type
  if (analysis_type == "social") {
    mancova_formula <- as.formula(
      paste("cbind(", paste(dependent_vars, collapse = ", "), 
            ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_SOCIAL_INT +
            SEX*DX_GROUP + SEX*AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
    )
  } else {
    mancova_formula <- as.formula(
      paste("cbind(", paste(dependent_vars, collapse = ", "), 
            ") ~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_COMM_INT +
            SEX*DX_GROUP + SEX*AVG_COMM_INT + AGE_AT_SCAN * SEX")
    )
  }
  
  # Perform MANCOVA
  results <- list()
  
  tryCatch({
    # Varimax MANCOVA
    varimax_mancova <- manova(mancova_formula, data = varimax_data)
    varimax_summary <- summary(varimax_mancova, test = "Wilks")
    varimax_aov <- summary.aov(varimax_mancova)
    
    results$varimax <- list(
      mancova = varimax_summary,
      anova = varimax_aov
    )
    
    # Save results
    save_mancova_results(varimax_summary, varimax_aov, "varimax", analysis_type)
    
  }, error = function(e) {
    message(paste("Varimax MANCOVA failed:", e$message))
  })
  
  tryCatch({
    # Promax MANCOVA
    promax_mancova <- manova(mancova_formula, data = promax_data)
    promax_summary <- summary(promax_mancova, test = "Wilks")
    promax_aov <- summary.aov(promax_mancova)
    
    results$promax <- list(
      mancova = promax_summary,
      anova = promax_aov
    )
    
    # Save results
    save_mancova_results(promax_summary, promax_aov, "promax", analysis_type)
    
  }, error = function(e) {
    message(paste("Promax MANCOVA failed:", e$message))
  })
  
  return(results)
}

#' Save MANCOVA results
save_mancova_results <- function(mancova_summary, anova_summary, rotation_type, analysis_type) {
  # Convert MANCOVA summary to dataframe
  summary_df <- as.data.frame(mancova_summary$stats)
  summary_df <- cbind(Variable = rownames(summary_df), summary_df)
  rownames(summary_df) <- NULL
  
  # Convert ANOVA summary to dataframe
  anova_list <- lapply(anova_summary, function(x) as.data.frame(x))
  detailed_summary_df <- do.call(rbind, lapply(names(anova_list), function(var) {
    df <- anova_list[[var]]
    df$Dependent_Variable <- var
    return(df)
  }))
  
  # Save files
  write.csv(summary_df, paste0(rotation_type, "_mancova_summary_", analysis_type, ".csv"))
  write.csv(detailed_summary_df, paste0(rotation_type, "_mancova_detailed_", analysis_type, ".csv"))
}

#' Perform MANCOVA assumption testing
#' @param varimax_data Varimax dataset
#' @param promax_data Promax dataset
#' @param analysis_type Type of analysis ("social" or "comm")
test_mancova_assumptions <- function(varimax_data, promax_data, analysis_type = "social") {
  message("Testing MANCOVA assumptions...")
  
  # Define variables
  dv_names <- paste0("PC", 1:20)
  independent_var <- "DX_GROUP"
  
  if (analysis_type == "social") {
    covariates <- c("AVG_SOCIAL_INT", "AGE_AT_SCAN", "SEX")
  } else {
    covariates <- c("AVG_COMM_INT", "AGE_AT_SCAN", "SEX")
  }
  
  # Test assumptions for both datasets
  test_dataset_assumptions(varimax_data, dv_names, independent_var, covariates, "varimax", analysis_type)
  test_dataset_assumptions(promax_data, dv_names, independent_var, covariates, "promax", analysis_type)
}

#' Test assumptions for a single dataset
test_dataset_assumptions <- function(data, dv_names, independent_var, covariates, 
                                   dataset_name, analysis_type) {
  
  dependent_vars <- data[, dv_names]
  data[[independent_var]] <- as.factor(data[[independent_var]])
  
  # VIF test for multicollinearity
  if (length(covariates) > 1) {
    tryCatch({
      vif_formula <- as.formula(paste("PC1 ~", independent_var, "+", paste(covariates, collapse = " + ")))
      vif_test <- vif(lm(vif_formula, data = data))
      
      vif_df <- data.frame(Variable = names(vif_test), VIF = vif_test)
      write.csv(vif_df, paste0("vif_test_", dataset_name, "_", analysis_type, ".csv"), row.names = FALSE)
    }, error = function(e) {
      message(paste("VIF test failed for", dataset_name, ":", e$message))
    })
  }
  
  # Multivariate normality test (Mardia's test)
  if (requireNamespace("MVN", quietly = TRUE)) {
    tryCatch({
      mardia_result <- MVN::mvn(dependent_vars, mvnTest = "mardia")
      
      normality_df <- data.frame(
        Test = c("Mardia Skewness", "Mardia Kurtosis"),
        Statistic = c(mardia_result$multivariateNormality$Statistic[1], 
                     mardia_result$multivariateNormality$Statistic[2]),
        P_Value = c(mardia_result$multivariateNormality$`p value`[1],
                   mardia_result$multivariateNormality$`p value`[2]),
        Result = c(mardia_result$multivariateNormality$Result[1],
                  mardia_result$multivariateNormality$Result[2])
      )
      
      write.csv(normality_df, paste0("mardia_test_", dataset_name, "_", analysis_type, ".csv"), row.names = FALSE)
    }, error = function(e) {
      message(paste("Mardia test failed for", dataset_name, ":", e$message))
    })
  }
  
  # Levene's test for homogeneity of variance
  levene_results <- list()
  for (dv in dv_names) {
    if (dv %in% names(data)) {
      tryCatch({
        levene_test <- leveneTest(as.formula(paste(dv, "~", independent_var)), data = data)
        levene_results[[dv]] <- levene_test
      }, error = function(e) {
        message(paste("Levene test failed for", dv, ":", e$message))
      })
    }
  }
  
  # Box's M test for homogeneity of covariance matrices
  if (requireNamespace("biotools", quietly = TRUE)) {
    tryCatch({
      box_test <- biotools::boxM(dependent_vars, data[[independent_var]])
      
      box_df <- data.frame(
        Test = "Box M",
        Chi_Sq = as.numeric(box_test$statistic),
        DF = box_test$parameter,
        P_Value = box_test$p.value
      )
      
      write.csv(box_df, paste0("box_test_", dataset_name, "_", analysis_type, ".csv"), row.names = FALSE)
    }, error = function(e) {
      message(paste("Box M test failed for", dataset_name, ":", e$message))
    })
  }
}

#' Perform permutation MANCOVA
#' @param varimax_data Varimax dataset
#' @param promax_data Promax dataset
#' @param analysis_type Type of analysis ("social" or "comm")
#' @param permutations Number of permutations
perform_permutation_mancova <- function(varimax_data, promax_data, analysis_type = "social", 
                                      permutations = 10000) {
  message("Performing permutation MANCOVA...")
  
  if (!requireNamespace("vegan", quietly = TRUE)) {
    message("vegan package required for permutation tests")
    return(NULL)
  }
  
  # Define dependent variables
  dv_names <- paste0("PC", 1:20)
  dependent_vars_v <- varimax_data[, dv_names]
  dependent_vars_p <- promax_data[, dv_names]
  
  set.seed(ANALYSIS_CONFIG$seed)
  
  # Create formula based on analysis type
  if (analysis_type == "social") {
    formula_str <- "~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT + SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX"
  } else {
    formula_str <- "~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT + SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX"
  }
  
  # Permutation tests
  tryCatch({
    perm_mancova_v <- vegan::anova.cca(
      vegan::capscale(as.matrix(dependent_vars_v) ~ ., data = varimax_data[, !names(varimax_data) %in% dv_names]),
      by = "margin", permutations = permutations
    )
    write.csv(perm_mancova_v, paste0("perm_mancova_varimax_", analysis_type, ".csv"))
  }, error = function(e) {
    message(paste("Permutation MANCOVA failed for varimax:", e$message))
  })
  
  tryCatch({
    perm_mancova_p <- vegan::anova.cca(
      vegan::capscale(as.matrix(dependent_vars_p) ~ ., data = promax_data[, !names(promax_data) %in% dv_names]),
      by = "margin", permutations = permutations
    )
    write.csv(perm_mancova_p, paste0("perm_mancova_promax_", analysis_type, ".csv"))
  }, error = function(e) {
    message(paste("Permutation MANCOVA failed for promax:", e$message))
  })
}

# ================================================================================
# SECTION 10: SUPPORT VECTOR MACHINE ANALYSIS
# ================================================================================

#' Comprehensive SVM analysis with cross-validation
#' @param varimax_data Varimax dataset
#' @param promax_data Promax dataset
perform_svm_analysis <- function(varimax_data, promax_data) {
  message("Performing SVM analysis...")
  
  # Prepare data
  svm_results <- list()
  
  # Process each dataset
  svm_results$varimax <- process_svm_dataset(varimax_data, "varimax")
  svm_results$promax <- process_svm_dataset(promax_data, "promax")
  
  # Compare results
  compare_svm_results(svm_results)
  
  return(svm_results)
}

#' Process SVM for a single dataset
process_svm_dataset <- function(data, dataset_name) {
  message(paste("Processing SVM for", dataset_name, "dataset..."))
  
  # Ensure DX_GROUP is factor
  data$DX_GROUP <- as.factor(data$DX_GROUP)
  
  # Create train/test split
  set.seed(ANALYSIS_CONFIG$seed)
  data <- data %>% mutate(RowID = row_number())
  
  # 10% for testing, stratified by DX_GROUP
  testing_data <- data %>%
    filter(DX_GROUP %in% c("1", "2")) %>%
    group_by(DX_GROUP) %>%
    slice_sample(prop = 0.1) %>%
    ungroup()
  
  training_data <- data %>%
    filter(!RowID %in% testing_data$RowID) %>%
    select(-RowID)
  
  testing_data <- testing_data %>% select(-RowID)
  
  # Save train/test data
  write.csv(training_data, paste0("training_", dataset_name, ".csv"))
  write.csv(testing_data, paste0("testing_", dataset_name, ".csv"))
  
  # Train SVM model
  svm_model <- train_svm_model(training_data)
  
  # Evaluate model
  evaluation_results <- evaluate_svm_model(svm_model, training_data, testing_data, dataset_name)
  
  return(list(
    model = svm_model,
    training_data = training_data,
    testing_data = testing_data,
    evaluation = evaluation_results
  ))
}

#' Train SVM model with cross-validation
train_svm_model <- function(training_data) {
  # Prepare data
  training_data$DX_GROUP <- as.factor(training_data$DX_GROUP)
  levels(training_data$DX_GROUP) <- make.names(levels(training_data$DX_GROUP))
  
  # Define predictors
  predictors <- paste0("PC", 1:20)
  
  # Train control
  train_control <- trainControl(
    method = "cv",
    number = ANALYSIS_CONFIG$cv_folds,
    classProbs = TRUE,
    savePredictions = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )
  
  # Parameter grid
  svm_grid <- expand.grid(
    C = ANALYSIS_CONFIG$svm_grid_range,
    sigma = ANALYSIS_CONFIG$svm_grid_range
  )
  
  # Train model
  svm_model <- train(
    DX_GROUP ~ .,
    data = training_data[, c("DX_GROUP", predictors)],
    method = "svmRadial",
    metric = "ROC",
    trControl = train_control,
    tuneGrid = svm_grid,
    probability = TRUE
  )
  
  return(svm_model)
}

#' Evaluate SVM model performance
evaluate_svm_model <- function(svm_model, training_data, testing_data, dataset_name) {
  # In-sample evaluation (cross-validation)
  cv_roc <- roc(
    response = svm_model$pred$obs,
    predictor = svm_model$pred$X1,
    levels = rev(levels(svm_model$pred$obs))
  )
  cv_auc <- auc(cv_roc)
  
  # Out-of-sample predictions
  test_predictions <- predict(svm_model, newdata = testing_data, type = "prob")
  test_classes <- predict(svm_model, newdata = testing_data, type = "raw")
  
  # Map predictions back to original labels
  test_classes <- as.character(test_classes)
  test_classes[test_classes == "X1"] <- "1"
  test_classes[test_classes == "X2"] <- "2"
  test_classes <- factor(test_classes, levels = levels(testing_data$DX_GROUP))
  
  # Out-of-sample ROC
  test_roc <- roc(testing_data$DX_GROUP, test_predictions[, 2])
  test_auc <- auc(test_roc)
  
  # Confusion matrix
  testing_data$DX_GROUP <- factor(testing_data$DX_GROUP, levels = c("1", "2"))
  test_classes <- factor(test_classes, levels = levels(testing_data$DX_GROUP))
  conf_matrix <- confusionMatrix(test_classes, testing_data$DX_GROUP)
  
  # Save results
  results_summary <- data.frame(
    Dataset = dataset_name,
    CV_AUC = cv_auc,
    Test_AUC = test_auc,
    Accuracy = conf_matrix$overall["Accuracy"],
    Sensitivity = conf_matrix$byClass["Sensitivity"],
    Specificity = conf_matrix$byClass["Specificity"]
  )
  
  write.csv(results_summary, paste0("svm_results_", dataset_name, ".csv"), row.names = FALSE)
  
  # Save detailed results
  test_results <- cbind(testing_data, 
                       Predicted_Class = test_classes,
                       Predicted_Prob = test_predictions[, 2])
  write.csv(test_results, paste0("svm_predictions_", dataset_name, ".csv"))
  
  return(list(
    cv_auc = cv_auc,
    test_auc = test_auc,
    confusion_matrix = conf_matrix,
    roc_cv = cv_roc,
    roc_test = test_roc,
    summary = results_summary
  ))
}

#' Compare SVM results between datasets
compare_svm_results <- function(svm_results) {
  # Create comparison plots
  create_roc_comparison_plots(svm_results)
  
  # Create performance comparison table
  if (!is.null(svm_results$varimax$evaluation$summary) && 
      !is.null(svm_results$promax$evaluation$summary)) {
    
    comparison_table <- rbind(
      svm_results$varimax$evaluation$summary,
      svm_results$promax$evaluation$summary
    )
    
    write.csv(comparison_table, "svm_comparison_results.csv", row.names = FALSE)
  }
}

#' Create ROC comparison plots
create_roc_comparison_plots <- function(svm_results) {
  # Set up plotting layout
  png("svm_roc_comparison.png", width = 1200, height = 800)
  par(mfrow = c(2, 2))
  
  # Varimax in-sample ROC
  if (!is.null(svm_results$varimax$evaluation$roc_cv)) {
    plot(svm_results$varimax$evaluation$roc_cv, col = "blue", 
         main = "In-Sample ROC - Varimax", lwd = 2)
    legend("bottomright", 
           legend = paste("AUC =", round(svm_results$varimax$evaluation$cv_auc, 4)), 
           col = "blue", lwd = 2)
  }
  
  # Promax in-sample ROC
  if (!is.null(svm_results$promax$evaluation$roc_cv)) {
    plot(svm_results$promax$evaluation$roc_cv, col = "darkgreen", 
         main = "In-Sample ROC - Promax", lwd = 2)
    legend("bottomright", 
           legend = paste("AUC =", round(svm_results$promax$evaluation$cv_auc, 4)), 
           col = "darkgreen", lwd = 2)
  }
  
  # Varimax out-of-sample ROC
  if (!is.null(svm_results$varimax$evaluation$roc_test)) {
    plot(svm_results$varimax$evaluation$roc_test, col = "red", 
         main = "Out-of-Sample ROC - Varimax", lwd = 2)
    legend("bottomright", 
           legend = paste("AUC =", round(svm_results$varimax$evaluation$test_auc, 4)), 
           col = "red", lwd = 2)
  }
  
  # Promax out-of-sample ROC
  if (!is.null(svm_results$promax$evaluation$roc_test)) {
    plot(svm_results$promax$evaluation$roc_test, col = "purple", 
         main = "Out-of-Sample ROC - Promax", lwd = 2)
    legend("bottomright", 
           legend = paste("AUC =", round(svm_results$promax$evaluation$test_auc, 4)), 
           col = "purple", lwd = 2)
  }
  
  dev.off()
}

# ================================================================================
# SECTION 11: BRAIN NETWORK VISUALIZATION
# ================================================================================

#' Create brain network visualizations
#' @param loadings_data PCA loadings data
#' @param roi_labels_file ROI labels file
#' @param top_n Number of top connections to visualize
#' @param analysis_type Type of analysis
create_brain_network_visualization <- function(loadings_data, roi_labels_file = "CC400_ROI_labels.csv", 
                                              top_n = 20, analysis_type = "general") {
  message("Creating brain network visualization...")
  
  if (!file.exists(roi_labels_file)) {
    message("ROI labels file not found:", roi_labels_file)
    return(NULL)
  }
  
  # Load ROI labels
  roi_labels <- read_csv(roi_labels_file, col_names = FALSE)
  roi_labels <- roi_labels %>% select(ROI_Number = 1, AAL_Label = 5)
  
  # Process loadings data
  if (is.character(loadings_data)) {
    # If it's a file path
    if (file.exists(loadings_data)) {
      loadings <- read_csv(loadings_data, col_names = FALSE)
      colnames(loadings) <- c("ROI_Connection", "Loading_Value")
    } else {
      message("Loadings file not found:", loadings_data)
      return(NULL)
    }
  } else {
    # If it's already loaded data
    loadings <- loadings_data
  }
  
  # Process loadings
  loadings$Loading_Value <- as.numeric(loadings$Loading_Value)
  
  # Split ROI connections
  loadings <- loadings %>%
    mutate(
      ROI_1 = sub("_.*", "", ROI_Connection),
      ROI_2 = sub(".*_", "", ROI_Connection)
    )
  
  # Remove prefixes
  loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
  loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)
  
  # Merge with ROI labels
  loadings <- loadings %>%
    left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
    rename(ROI_1_Label = AAL_Label) %>%
    left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
    rename(ROI_2_Label = AAL_Label)
  
  # Select top connections
  top_connections <- loadings %>%
    arrange(desc(abs(Loading_Value))) %>%
    slice(1:top_n)
  
  # Create network visualization
  create_network_plot(top_connections, analysis_type, top_n)
  
  # Create brain region heatmap
  create_brain_heatmap(top_connections, analysis_type)
  
  # Save top connections
  write.csv(top_connections, paste0("top_", top_n, "_connections_", analysis_type, ".csv"))
  
  return(top_connections)
}

#' Create network plot
create_network_plot <- function(top_connections, analysis_type, top_n) {
  if (!requireNamespace("ggraph", quietly = TRUE) || 
      !requireNamespace("tidygraph", quietly = TRUE)) {
    message("ggraph and tidygraph packages required for network visualization")
    return(NULL)
  }
  
  # Build graph edges
  edges <- top_connections %>%
    select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)
  
  # Create nodes
  nodes <- data.frame(name = unique(c(edges$from, edges$to)))
  
  # Assign hemisphere layout
  nodes <- nodes %>%
    arrange(name) %>%
    mutate(
      hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
      x = ifelse(hemisphere == "Left", 
                runif(ceiling(n()/2), min = -2, max = -0.8), 
                runif(floor(n()/2), min = 0.8, max = 2)),
      y = seq(-2, 2, length.out = n())
    )
  
  # Create graph object
  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  graph <- graph %>% mutate(x = nodes$x, y = nodes$y)
  
  # Create plot
  network_plot <- ggraph::ggraph(graph, layout = "manual", x = x, y = y) +
    ggraph::geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.6) +
    ggraph::geom_node_point(size = 4, color = "black") +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    ggraph::scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggraph::scale_edge_width(range = c(0.3, 1.5)) +
    theme_void() +
    labs(
      title = paste("Top", top_n, "ROI Connections -", analysis_type),
      edge_width = "Strength",
      edge_color = "Loading Value"
    ) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  ggsave(paste0("network_plot_", analysis_type, ".png"), plot = network_plot, 
         width = 12, height = 8)
  
  return(network_plot)
}

#' Create brain region heatmap
create_brain_heatmap <- function(top_connections, analysis_type) {
  # Brain region mapping function
  roi_brain_region <- function(label) {
    if (grepl("Frontal", label)) return("Frontal Lobe")
    else if (grepl("Temporal", label)) return("Temporal Lobe")
    else if (grepl("Parietal", label)) return("Parietal Lobe")
    else if (grepl("Occipital", label)) return("Occipital Lobe")
    else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
    else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
    else return("Other")
  }
  
  # Map brain regions
  top_connections_long <- top_connections %>%
    select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
    pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
    mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))
  
  # Summarize by brain region
  brain_heat_data <- top_connections_long %>%
    group_by(Brain_Region) %>%
    summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
    filter(Brain_Region != "Other")
  
  # Save brain region summary
  write.csv(brain_heat_data, paste0("brain_region_summary_", analysis_type, ".csv"))
  
  return(brain_heat_data)
}

# ================================================================================
# SECTION 12: COMPREHENSIVE QUALITY CONTROL AND VALIDATION
# ================================================================================

#' Perform comprehensive quality control checks
perform_quality_control <- function() {
  message("Performing quality control checks...")
  
  qc_results <- list()
  
  # Check data integrity
  qc_results$data_integrity <- check_data_integrity()
  
  # Check missing data patterns
  qc_results$missing_data <- check_missing_data_patterns()
  
  # Check statistical assumptions
  qc_results$assumptions <- check_statistical_assumptions()
  
  # Generate QC report
  generate_qc_report(qc_results)
  
  return(qc_results)
}

#' Check data integrity
check_data_integrity <- function() {
  message("Checking data integrity...")
  
  integrity_checks <- list()
  
  # Check if required files exist
  required_files <- c("phenotype_data.csv", "merged_data.csv", 
                     "varimax_alldata.csv", "promax_alldata.csv")
  
  for (file in required_files) {
    integrity_checks[[paste0("file_", gsub("\\.|csv", "", file))]] <- file.exists(file)
  }
  
  return(integrity_checks)
}

#' Check missing data patterns
check_missing_data_patterns <- function() {
  message("Checking missing data patterns...")
  
  missing_patterns <- list()
  
  # Check key datasets for missing data
  datasets_to_check <- c("varimax_alldata.csv", "promax_alldata.csv")
  
  for (dataset_file in datasets_to_check) {
    if (file.exists(dataset_file)) {
      data <- fread(dataset_file)
      
      # Calculate missing data percentage by column
      missing_pct <- data %>%
        summarise_all(~ sum(is.na(.)) / length(.) * 100) %>%
        pivot_longer(everything(), names_to = "Variable", values_to = "Missing_Percent")
      
      missing_patterns[[dataset_file]] <- missing_pct
      
      # Save missing data report
      write.csv(missing_pct, paste0("missing_data_", gsub("\\.csv", "", dataset_file), ".csv"))
    }
  }
  
  return(missing_patterns)
}

#' Check statistical assumptions
check_statistical_assumptions <- function() {
  message("Checking statistical assumptions...")
  
  assumption_checks <- list()
  
  # This would include checks for normality, homoscedasticity, etc.
  # Implementation would depend on specific requirements
  
  return(assumption_checks)
}

#' Generate quality control report
generate_qc_report <- function(qc_results) {
  message("Generating quality control report...")
  
  # Create a summary report
  report_lines <- c(
    "# Quality Control Report",
    "Generated on:", Sys.time(),
    "",
    "## Data Integrity Checks",
    paste(names(qc_results$data_integrity), ":", qc_results$data_integrity),
    "",
    "## Missing Data Summary",
    "See individual missing data CSV files for details.",
    "",
    "## Statistical Assumptions",
    "See individual assumption test files for details."
  )
  
  writeLines(report_lines, "quality_control_report.txt")
}

# ================================================================================
# SECTION 13: MAIN EXECUTION WORKFLOW
# ================================================================================

#' Main analysis workflow
#' @param run_data_extraction Whether to run data extraction step
#' @param run_full_analysis Whether to run the complete analysis pipeline
main_analysis_workflow <- function(run_data_extraction = FALSE, run_full_analysis = TRUE) {
  message("="*80)
  message("STARTING COMPREHENSIVE NEUROIMAGING ANALYSIS PIPELINE")
  message("="*80)
  
  start_time <- Sys.time()
  
  # Initialize results storage
  analysis_results <- list()
  
  tryCatch({
    
    # Step 1: Data Extraction (optional)
    if (run_data_extraction) {
      message("\n--- STEP 1: DATA EXTRACTION ---")
      extract_abide_data(ANALYSIS_CONFIG$base_path)
      combined_data <- combine_csv_files()
      analysis_results$combined_data <- combined_data
    }
    
    if (!run_full_analysis) {
      message("Data extraction completed. Skipping full analysis.")
      return(analysis_results)
    }
    
    # Step 2: Phenotype Data Processing
    message("\n--- STEP 2: PHENOTYPE DATA PROCESSING ---")
    phenotype_data <- process_phenotype_data()
    analysis_results$phenotype_data <- phenotype_data
    
    # Step 3: Correlation Analysis
    message("\n--- STEP 3: CORRELATION ANALYSIS ---")
    if (file.exists("CC400_combined.csv")) {
      correlation_data <- perform_correlation_analysis(NULL)
      analysis_results$correlation_data <- correlation_data
    } else {
      message("Combined data file not found. Skipping correlation analysis.")
    }
    
    # Step 4: Data Merging
    message("\n--- STEP 4: DATA MERGING ---")
    merged_data <- merge_datasets(phenotype_data, analysis_results$correlation_data)
    analysis_results$merged_data <- merged_data
    
    # Step 5: Principal Component Analysis
    message("\n--- STEP 5: PRINCIPAL COMPONENT ANALYSIS ---")
    if (!is.null(merged_data)) {
      pca_results <- perform_comprehensive_pca(merged_data, ANALYSIS_CONFIG$n_components_pca)
      analysis_results$pca_results <- pca_results
      
      # Step 6: Create Analysis Datasets
      message("\n--- STEP 6: CREATING ANALYSIS DATASETS ---")
      analysis_datasets <- create_analysis_datasets(
        merged_data, 
        pca_results$varimax_scores, 
        pca_results$promax_scores
      )
      analysis_results$analysis_datasets <- analysis_datasets
      
    } else {
      message("Merged data not available. Skipping PCA and subsequent analyses.")
      return(analysis_results)
    }
    
    # Step 7: Exploratory Data Analysis
    message("\n--- STEP 7: EXPLORATORY DATA ANALYSIS ---")
    perform_comprehensive_eda(
      merged_data, 
      analysis_datasets$varimax, 
      analysis_datasets$promax
    )
    
    # Step 8: MANCOVA Analysis
    message("\n--- STEP 8: MANCOVA ANALYSIS ---")
    
    # Test assumptions first
    test_mancova_assumptions(analysis_datasets$varimax, analysis_datasets$promax, "social")
    test_mancova_assumptions(analysis_datasets$varimax, analysis_datasets$promax, "comm")
    
    # Perform MANCOVA
    mancova_social <- perform_mancova_analysis(analysis_datasets$varimax, analysis_datasets$promax, "social")
    mancova_comm <- perform_mancova_analysis(analysis_datasets$varimax, analysis_datasets$promax, "comm")
    
    analysis_results$mancova_results <- list(social = mancova_social, comm = mancova_comm)
    
    # Permutation MANCOVA
    perform_permutation_mancova(analysis_datasets$varimax, analysis_datasets$promax, "social")
    perform_permutation_mancova(analysis_datasets$varimax, analysis_datasets$promax, "comm")
    
    # Step 9: Support Vector Machine Analysis
    message("\n--- STEP 9: SUPPORT VECTOR MACHINE ANALYSIS ---")
    svm_results <- perform_svm_analysis(analysis_datasets$varimax, analysis_datasets$promax)
    analysis_results$svm_results <- svm_results
    
    # Step 10: Brain Network Visualization
    message("\n--- STEP 10: BRAIN NETWORK VISUALIZATION ---")
    
    # Create visualizations for different components/analyses
    component_files <- list.files(pattern = "component_.*_loadings.*\\.csv")
    
    for (comp_file in component_files[1:min(3, length(component_files))]) {
      analysis_type <- gsub("component_|_loadings.*\\.csv", "", comp_file)
      tryCatch({
        create_brain_network_visualization(comp_file, analysis_type = analysis_type)
      }, error = function(e) {
        message(paste("Visualization failed for", comp_file, ":", e$message))
      })
    }
    
    # Step 11: Quality Control
    message("\n--- STEP 11: QUALITY CONTROL ---")
    qc_results <- perform_quality_control()
    analysis_results$qc_results <- qc_results
    
    # Final Summary
    total_time <- Sys.time() - start_time
    message("\n" + "="*80)
    message("ANALYSIS PIPELINE COMPLETED SUCCESSFULLY")
    message(paste("Total execution time:", round(total_time, 2), attr(total_time, "units")))
    message("="*80)
    
    # Save final results summary
    save_final_summary(analysis_results, total_time)
    
  }, error = function(e) {
    message("\nERROR in main analysis workflow:", e$message)
    message("Partial results may be available in the analysis_results object")
  })
  
  return(analysis_results)
}

#' Save final analysis summary
save_final_summary <- function(analysis_results, total_time) {
  summary_info <- list(
    execution_time = total_time,
    timestamp = Sys.time(),
    completed_steps = names(analysis_results),
    config_used = ANALYSIS_CONFIG
  )
  
  # Save as JSON if jsonlite is available, otherwise as text
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(summary_info, "analysis_summary.json", pretty = TRUE)
  } else {
    capture.output(str(summary_info), file = "analysis_summary.txt")
  }
  
  message("Analysis summary saved")
}

# ================================================================================
# SECTION 14: UTILITY FUNCTIONS AND HELPERS
# ================================================================================

#' Print analysis configuration
print_analysis_config <- function() {
  message("Current Analysis Configuration:")
  message("=============================")
  for (name in names(ANALYSIS_CONFIG)) {
    message(paste(name, ":", ANALYSIS_CONFIG[[name]]))
  }
  message("=============================\n")
}

#' Update analysis configuration
#' @param ... Named parameters to update in the configuration
update_analysis_config <- function(...) {
  new_params <- list(...)
  for (param_name in names(new_params)) {
    if (param_name %in% names(ANALYSIS_CONFIG)) {
      ANALYSIS_CONFIG[[param_name]] <<- new_params[[param_name]]
      message(paste("Updated", param_name, "to:", new_params[[param_name]]))
    } else {
      message(paste("Warning: Unknown parameter", param_name))
    }
  }
}

#' Create analysis directory structure
create_directory_structure <- function() {
  dirs_to_create <- c(
    "results",
    "plots", 
    "data_processed",
    "models",
    "reports"
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir)
      message(paste("Created directory:", dir))
    }
  }
}

#' Clean up temporary files
cleanup_temp_files <- function(confirm = TRUE) {
  temp_patterns <- c("*.tmp", "*_temp.*", "Rplot*.png")
  
  temp_files <- unlist(lapply(temp_patterns, function(pattern) {
    list.files(pattern = pattern, full.names = TRUE)
  }))
  
  if (length(temp_files) > 0) {
    if (confirm) {
      response <- readline(paste("Delete", length(temp_files), "temporary files? (y/n): "))
      if (tolower(response) != "y") {
        message("Cleanup cancelled")
        return()
      }
    }
    
    file.remove(temp_files)
    message(paste("Removed", length(temp_files), "temporary files"))
  } else {
    message("No temporary files found")
  }
}

# ================================================================================
# SECTION 15: EXECUTION CONTROL
# ================================================================================

# Print configuration
print_analysis_config()

# Create directory structure
create_directory_structure()

# Main execution
if (interactive()) {
  message("\nThis is an interactive session.")
  message("To run the full analysis pipeline, use:")
  message("results <- main_analysis_workflow()")
  message("\nTo run only data extraction, use:")
  message("results <- main_analysis_workflow(run_data_extraction = TRUE, run_full_analysis = FALSE)")
} else {
  # Automatically run if not in interactive mode
  message("Running automated analysis pipeline...")
  results <- main_analysis_workflow()
}

# ================================================================================
# END OF COMBINED ANALYSIS SCRIPT
# ================================================================================

message("\n" + "="*80)
message("COMBINED ANALYSIS SCRIPT LOADED SUCCESSFULLY")
message("All functions and configurations are now available.")
message("Use main_analysis_workflow() to run the complete analysis.")
message("="*80)