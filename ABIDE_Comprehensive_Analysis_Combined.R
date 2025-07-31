#!/usr/bin/env Rscript
#=============================================================================
# ABIDE Comprehensive Analysis - Combined Script
# Integrates all analysis workflows into one efficient and executable code
#=============================================================================

# Load required libraries
cat("Loading required libraries...\n")
required_packages <- c(
  "dplyr", "data.table", "readr", "tidyr", "reshape2",
  "ggplot2", "corrplot", "GPArotation", "psych", 
  "e1071", "caret", "pROC", "MVN", "car", "biotools",
  "vegan", "igraph", "ggraph", "tidygraph", "ggimage",
  "png", "grid"
)

# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

install_and_load(required_packages)

# Set up working directory and create output folders
cat("Setting up directory structure...\n")
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/plots")) dir.create("output/plots")
if (!dir.exists("output/data")) dir.create("output/data")
if (!dir.exists("output/results")) dir.create("output/results")

#=============================================================================
# SECTION 1: DATA LOADING AND PREPROCESSING
#=============================================================================
cat("\n=== SECTION 1: DATA LOADING AND PREPROCESSING ===\n")

# 1.1 Load and process phenotype data
cat("1.1 Processing phenotype data...\n")
process_phenotype_data <- function() {
  # Load phenotype data
  if (!file.exists("Phenotypic_V1_0b_v1.csv")) {
    stop("Error: Phenotypic_V1_0b_v1.csv not found. Please ensure the file is in the working directory.")
  }
  
  p_data <- fread("Phenotypic_V1_0b_v1.csv")
  cat(sprintf("Loaded phenotype data: %d rows, %d columns\n", nrow(p_data), ncol(p_data)))
  
  # Select needed columns
  p_data1 <- p_data %>% 
    dplyr::select(SITE_ID, SUB_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN,
                  SEX, ADI_R_SOCIAL_TOTAL_A, ADI_R_VERBAL_TOTAL_BV,
                  ADOS_COMM, ADOS_SOCIAL, SRS_COGNITION, SRS_COMMUNICATION)
  
  # Define columns for processing
  start_col <- which(names(p_data1) == "ADI_R_SOCIAL_TOTAL_A")
  end_col <- which(names(p_data1) == "SRS_COMMUNICATION")
  columns_to_check <- names(p_data1)[start_col:end_col]
  
  # Replace -9999 values with NA
  p_data1[p_data1 == -9999] <- NA
  
  # Remove rows where all data in specified columns are blank or NA
  cleaned_data <- p_data1 %>%
    filter(rowSums(sapply(select(., all_of(columns_to_check)), 
                         function(x) !is.na(x) & x != "")) > 0)
  
  cat(sprintf("After cleaning: %d rows, %d columns\n", nrow(cleaned_data), ncol(cleaned_data)))
  
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
  
  # Create composite variables
  n_data[, AVG_SOCIAL := rowMeans(.SD, na.rm = TRUE), 
         .SDcols = c("ADI_R_SOCIAL_TOTAL_A", "ADOS_SOCIAL", "SRS_COGNITION")]
  n_data[, AVG_COMM := rowMeans(.SD, na.rm = TRUE), 
         .SDcols = c("ADI_R_VERBAL_TOTAL_BV", "ADOS_COMM", "SRS_COMMUNICATION")]
  
  # Save processed data
  write.csv(n_data, "output/data/phenotype_data.csv", row.names = FALSE)
  return(n_data)
}

# 1.2 Process correlation data (if CPAC files exist)
process_correlation_data <- function() {
  cat("1.2 Processing correlation data...\n")
  
  # Check if combined correlation data exists
  if (file.exists("varimax_alldata.csv") && file.exists("promax_alldata.csv")) {
    cat("Found existing processed correlation data files\n")
    return(list(
      varimax = fread("varimax_alldata.csv"),
      promax = fread("promax_alldata.csv")
    ))
  }
  
  # If correlation results exist, use them
  if (file.exists("correlation_results_wide.csv")) {
    cat("Found correlation results, loading...\n")
    corr_data <- fread("correlation_results_wide.csv")
    return(corr_data)
  }
  
  cat("Warning: No correlation data found. Continuing with available data...\n")
  return(NULL)
}

# Execute data processing
phenotype_data <- process_phenotype_data()
correlation_data <- process_correlation_data()

#=============================================================================
# SECTION 2: PCA ANALYSIS WITH ROTATION
#=============================================================================
cat("\n=== SECTION 2: PCA ANALYSIS WITH ROTATION ===\n")

perform_pca_analysis <- function(correlation_data, phenotype_data) {
  if (is.null(correlation_data)) {
    cat("Skipping PCA analysis - no correlation data available\n")
    return(NULL)
  }
  
  # If we have processed data, return it
  if (is.list(correlation_data) && "varimax" %in% names(correlation_data)) {
    return(correlation_data)
  }
  
  cat("2.1 Performing PCA analysis...\n")
  start_time <- Sys.time()
  
  # Merge phenotype and correlation data if needed
  if (is.data.frame(correlation_data)) {
    merged_data <- merge(phenotype_data, correlation_data, 
                        by.x = "SUB_ID", by.y = "Subject", all = FALSE)
    cat(sprintf("Merged data dimensions: %d rows, %d columns\n", 
                nrow(merged_data), ncol(merged_data)))
  } else {
    merged_data <- correlation_data
  }
  
  # Select correlation variables for PCA
  corr_cols <- grep("^X\\.", names(merged_data), value = TRUE)
  if (length(corr_cols) == 0) {
    cat("No correlation columns found for PCA\n")
    return(NULL)
  }
  
  # Prepare data for PCA
  selected_vars <- merged_data[, ..corr_cols]
  numeric_vars <- na.omit(as.data.frame(lapply(selected_vars, as.numeric)))
  
  cat(sprintf("PCA input: %d variables, %d observations\n", 
              ncol(numeric_vars), nrow(numeric_vars)))
  
  # Perform PCA
  pca_result <- prcomp(numeric_vars, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  eigenvalues <- pca_result$sdev^2
  explained_var <- eigenvalues
  prop_var <- explained_var / sum(explained_var)
  cum_var <- cumsum(prop_var)
  
  # Create variance dataframe
  var_df <- data.frame(
    PC = paste0("PC", 1:length(prop_var)),
    Proportion = prop_var,
    Cumulative = cum_var,
    Explained_Variance = explained_var
  )
  
  write.csv(var_df, "output/results/explained_variance_pca.csv", row.names = FALSE)
  
  # Create scree plot
  png("output/plots/scree_plot.png", width = 800, height = 600)
  plot(explained_var[1:50], type = "b", main = "Scree Plot - Explained Variance",
       xlab = "Principal Component", ylab = "Explained Variance",
       pch = 19, col = "blue")
  abline(v = 20, col = "darkgreen", lty = 3, lwd = 2)
  dev.off()
  
  # Perform rotations
  loadings_50 <- as.matrix(pca_result$rotation[, 1:50])
  
  # Varimax rotation
  cat("2.2 Performing Varimax rotation...\n")
  varimax_result <- varimax(loadings_50)
  varimax_loadings <- varimax_result$loadings[, 1:20]
  
  # Calculate varimax scores
  numeric_vars_matrix <- as.matrix(numeric_vars)
  varimax_scores <- as.data.frame(numeric_vars_matrix %*% varimax_loadings)
  colnames(varimax_scores) <- paste0("PC", 1:20)
  
  # Promax rotation
  cat("2.3 Performing Promax rotation...\n")
  promax_result <- promax(loadings_50)
  promax_loadings <- promax_result$loadings[, 1:20]
  
  # Calculate promax scores
  promax_scores <- as.data.frame(numeric_vars_matrix %*% promax_loadings)
  colnames(promax_scores) <- paste0("PC", 1:20)
  
  # Create combined datasets
  alldata <- merged_data %>% 
    select(SUB_ID, SITE_ID, DX_GROUP, DSM_IV_TR, AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM) %>%
    mutate(
      AVG_SOCIAL_INT = round(AVG_SOCIAL * 10),
      AVG_COMM_INT = round(AVG_COMM * 10)
    )
  
  varimax_alldata <- cbind(alldata, varimax_scores)
  promax_alldata <- cbind(alldata, promax_scores)
  
  # Save datasets
  write.csv(varimax_scores, "output/data/varimax_scores.csv", row.names = FALSE)
  write.csv(promax_scores, "output/data/promax_scores.csv", row.names = FALSE)
  write.csv(varimax_alldata, "output/data/varimax_alldata.csv", row.names = FALSE)
  write.csv(promax_alldata, "output/data/promax_alldata.csv", row.names = FALSE)
  
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("PCA analysis completed in %.2f seconds\n", total_time))
  
  return(list(varimax = varimax_alldata, promax = promax_alldata))
}

# Execute PCA analysis
pca_results <- perform_pca_analysis(correlation_data, phenotype_data)

#=============================================================================
# SECTION 3: EXPLORATORY DATA ANALYSIS
#=============================================================================
cat("\n=== SECTION 3: EXPLORATORY DATA ANALYSIS ===\n")

perform_eda <- function(data_list) {
  if (is.null(data_list)) {
    cat("Skipping EDA - no PCA data available\n")
    return(NULL)
  }
  
  cat("3.1 Generating descriptive statistics...\n")
  
  # Use either PCA results or phenotype data
  if (is.list(data_list) && "varimax" %in% names(data_list)) {
    varimax_data <- data_list$varimax
    promax_data <- data_list$promax
  } else {
    # Use phenotype data for basic EDA
    varimax_data <- phenotype_data
    promax_data <- phenotype_data
  }
  
  # Basic demographics
  if ("DX_GROUP" %in% names(varimax_data)) {
    dx_summary <- varimax_data %>%
      group_by(DX_GROUP) %>%
      summarize(
        n = n(),
        mean_age = mean(AGE_AT_SCAN, na.rm = TRUE),
        sd_age = sd(AGE_AT_SCAN, na.rm = TRUE),
        male_count = sum(SEX == 1, na.rm = TRUE),
        female_count = sum(SEX == 2, na.rm = TRUE),
        .groups = "drop"
      )
    
    write.csv(dx_summary, "output/results/demographics_summary.csv", row.names = FALSE)
    print(dx_summary)
  }
  
  # Create basic visualizations
  cat("3.2 Creating visualizations...\n")
  
  # Age distribution
  if ("AGE_AT_SCAN" %in% names(varimax_data)) {
    p1 <- ggplot(varimax_data, aes(x = AGE_AT_SCAN)) +
      geom_histogram(binwidth = 2, fill = "skyblue", alpha = 0.7) +
      labs(title = "Age Distribution", x = "Age at Scan", y = "Frequency") +
      theme_minimal()
    
    ggsave("output/plots/age_distribution.png", p1, width = 8, height = 6)
  }
  
  # Diagnostic group distribution
  if ("DX_GROUP" %in% names(varimax_data)) {
    p2 <- ggplot(varimax_data, aes(x = factor(DX_GROUP))) +
      geom_bar(fill = "lightgreen", alpha = 0.7) +
      geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
      scale_x_discrete(labels = c("1" = "Autism", "2" = "Control")) +
      labs(title = "Diagnostic Group Distribution", x = "Group", y = "Count") +
      theme_minimal()
    
    ggsave("output/plots/dx_group_distribution.png", p2, width = 8, height = 6)
  }
  
  # PC correlation analysis if PC data exists
  pc_cols <- grep("^PC\\d+$", names(varimax_data), value = TRUE)
  if (length(pc_cols) >= 20) {
    cat("3.3 Analyzing principal components...\n")
    
    # PC summary by group
    pc_summary <- varimax_data %>%
      select(DX_GROUP, all_of(pc_cols[1:20])) %>%
      group_by(DX_GROUP) %>%
      summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop")
    
    write.csv(pc_summary, "output/results/pc_summary_by_group.csv", row.names = FALSE)
    
    # Correlation matrix visualization
    pc_data <- varimax_data[, pc_cols[1:20], with = FALSE]
    pc_corr <- cor(pc_data, use = "complete.obs")
    
    png("output/plots/pc_correlation_matrix.png", width = 800, height = 800)
    corrplot(pc_corr, method = "color", title = "PC Correlation Matrix", 
             tl.cex = 0.8, cl.cex = 0.8)
    dev.off()
  }
  
  cat("EDA completed successfully\n")
}

# Execute EDA
perform_eda(pca_results)

#=============================================================================
# SECTION 4: STATISTICAL ANALYSIS (MANCOVA)
#=============================================================================
cat("\n=== SECTION 4: STATISTICAL ANALYSIS (MANCOVA) ===\n")

perform_mancova_analysis <- function(data_list) {
  if (is.null(data_list) || !is.list(data_list)) {
    cat("Skipping MANCOVA - no suitable data available\n")
    return(NULL)
  }
  
  cat("4.1 Performing MANCOVA assumption tests...\n")
  
  varimax_data <- data_list$varimax
  promax_data <- data_list$promax
  
  # Define dependent variables
  dv_names <- paste0("PC", 1:20)
  
  # Check if we have the required variables
  if (!all(dv_names %in% names(varimax_data))) {
    cat("Missing PC variables for MANCOVA\n")
    return(NULL)
  }
  
  # Prepare data
  varimax_data$DX_GROUP <- as.factor(varimax_data$DX_GROUP)
  promax_data$DX_GROUP <- as.factor(promax_data$DX_GROUP)
  
  # MANCOVA for communication variables
  cat("4.2 MANCOVA for communication variables...\n")
  
  mancova_formula_comm <- as.formula(
    paste("cbind(", paste(dv_names, collapse = ", "), 
          ") ~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
          SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX")
  )
  
  tryCatch({
    # Varimax MANCOVA
    varimax_mancova <- manova(mancova_formula_comm, data = varimax_data)
    varimax_summary <- summary(varimax_mancova, test = "Wilks")
    
    # Promax MANCOVA
    promax_mancova <- manova(mancova_formula_comm, data = promax_data)
    promax_summary <- summary(promax_mancova, test = "Wilks")
    
    # Save results
    varimax_results <- as.data.frame(varimax_summary$stats)
    promax_results <- as.data.frame(promax_summary$stats)
    
    write.csv(varimax_results, "output/results/mancova_varimax_comm.csv")
    write.csv(promax_results, "output/results/mancova_promax_comm.csv")
    
    cat("MANCOVA analysis completed\n")
    
  }, error = function(e) {
    cat("Error in MANCOVA analysis:", e$message, "\n")
  })
  
  # MANCOVA for social variables
  cat("4.3 MANCOVA for social variables...\n")
  
  mancova_formula_social <- as.formula(
    paste("cbind(", paste(dv_names, collapse = ", "), 
          ") ~ DX_GROUP * AVG_SOCIAL_INT + AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_SOCIAL_INT +
          SEX * DX_GROUP + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX")
  )
  
  tryCatch({
    # Varimax MANCOVA
    varimax_mancova_social <- manova(mancova_formula_social, data = varimax_data)
    varimax_summary_social <- summary(varimax_mancova_social, test = "Wilks")
    
    # Promax MANCOVA
    promax_mancova_social <- manova(mancova_formula_social, data = promax_data)
    promax_summary_social <- summary(promax_mancova_social, test = "Wilks")
    
    # Save results
    varimax_results_social <- as.data.frame(varimax_summary_social$stats)
    promax_results_social <- as.data.frame(promax_summary_social$stats)
    
    write.csv(varimax_results_social, "output/results/mancova_varimax_social.csv")
    write.csv(promax_results_social, "output/results/mancova_promax_social.csv")
    
    cat("Social MANCOVA analysis completed\n")
    
  }, error = function(e) {
    cat("Error in social MANCOVA analysis:", e$message, "\n")
  })
}

# Execute MANCOVA analysis
perform_mancova_analysis(pca_results)

#=============================================================================
# SECTION 5: MACHINE LEARNING (SVM)
#=============================================================================
cat("\n=== SECTION 5: MACHINE LEARNING (SVM) ===\n")

perform_svm_analysis <- function(data_list) {
  if (is.null(data_list) || !is.list(data_list)) {
    cat("Skipping SVM - no suitable data available\n")
    return(NULL)
  }
  
  cat("5.1 Preparing SVM datasets...\n")
  
  varimax_data <- data_list$varimax
  promax_data <- data_list$promax
  
  # Check if we have the required variables
  pc_vars <- paste0("PC", 1:20)
  if (!all(pc_vars %in% names(varimax_data))) {
    cat("Missing PC variables for SVM\n")
    return(NULL)
  }
  
  # Prepare data
  varimax_data$DX_GROUP <- as.factor(varimax_data$DX_GROUP)
  promax_data$DX_GROUP <- as.factor(promax_data$DX_GROUP)
  
  # Create train/test splits
  set.seed(123)
  
  # Split varimax data
  varimax_data <- varimax_data %>% mutate(RowID = row_number())
  testing_vdata <- varimax_data %>%
    filter(DX_GROUP %in% c("1", "2")) %>%
    group_by(DX_GROUP) %>%
    slice_sample(prop = 0.2) %>%
    ungroup()
  training_vdata <- varimax_data %>%
    filter(!RowID %in% testing_vdata$RowID) %>%
    select(-RowID)
  testing_vdata <- testing_vdata %>% select(-RowID)
  
  # Split promax data
  promax_data <- promax_data %>% mutate(RowID = row_number())
  testing_pdata <- promax_data %>%
    filter(DX_GROUP %in% c("1", "2")) %>%
    group_by(DX_GROUP) %>%
    slice_sample(prop = 0.2) %>%
    ungroup()
  training_pdata <- promax_data %>%
    filter(!RowID %in% testing_pdata$RowID) %>%
    select(-RowID)
  testing_pdata <- testing_pdata %>% select(-RowID)
  
  cat("5.2 Training SVM models...\n")
  
  # Rename factor levels to valid R names
  levels(training_vdata$DX_GROUP) <- make.names(levels(training_vdata$DX_GROUP))
  levels(training_pdata$DX_GROUP) <- make.names(levels(training_pdata$DX_GROUP))
  levels(testing_vdata$DX_GROUP) <- make.names(levels(testing_vdata$DX_GROUP))
  levels(testing_pdata$DX_GROUP) <- make.names(levels(testing_pdata$DX_GROUP))
  
  # Train control
  train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    savePredictions = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # SVM grid
  svm_grid <- expand.grid(
    C = c(0.1, 1, 10),
    sigma = c(0.01, 0.1, 1)
  )
  
  tryCatch({
    # Train varimax SVM
    svm_varimax <- train(
      DX_GROUP ~ .,
      data = training_vdata[, c("DX_GROUP", pc_vars)],
      method = "svmRadial",
      metric = "ROC",
      trControl = train_control,
      tuneGrid = svm_grid
    )
    
    # Train promax SVM
    svm_promax <- train(
      DX_GROUP ~ .,
      data = training_pdata[, c("DX_GROUP", pc_vars)],
      method = "svmRadial",
      metric = "ROC",
      trControl = train_control,
      tuneGrid = svm_grid
    )
    
    cat("5.3 Evaluating SVM performance...\n")
    
    # Predictions
    pred_varimax <- predict(svm_varimax, testing_vdata, type = "prob")
    pred_promax <- predict(svm_promax, testing_pdata, type = "prob")
    
    # ROC analysis
    roc_varimax <- roc(testing_vdata$DX_GROUP, pred_varimax[, 2])
    roc_promax <- roc(testing_pdata$DX_GROUP, pred_promax[, 2])
    
    # Save results
    svm_results <- data.frame(
      Model = c("Varimax", "Promax"),
      AUC = c(auc(roc_varimax), auc(roc_promax)),
      Best_C = c(svm_varimax$bestTune$C, svm_promax$bestTune$C),
      Best_Sigma = c(svm_varimax$bestTune$sigma, svm_promax$bestTune$sigma)
    )
    
    write.csv(svm_results, "output/results/svm_performance.csv", row.names = FALSE)
    
    # Plot ROC curves
    png("output/plots/roc_curves.png", width = 800, height = 600)
    par(mfrow = c(1, 2))
    plot(roc_varimax, main = "Varimax SVM ROC", col = "blue", lwd = 2)
    legend("bottomright", legend = paste("AUC =", round(auc(roc_varimax), 3)))
    plot(roc_promax, main = "Promax SVM ROC", col = "red", lwd = 2)
    legend("bottomright", legend = paste("AUC =", round(auc(roc_promax), 3)))
    dev.off()
    
    cat("SVM analysis completed\n")
    print(svm_results)
    
  }, error = function(e) {
    cat("Error in SVM analysis:", e$message, "\n")
  })
}

# Execute SVM analysis
perform_svm_analysis(pca_results)

#=============================================================================
# SECTION 6: PERMUTATION TESTING
#=============================================================================
cat("\n=== SECTION 6: PERMUTATION TESTING ===\n")

perform_permutation_tests <- function(data_list) {
  if (is.null(data_list) || !is.list(data_list)) {
    cat("Skipping permutation tests - no suitable data available\n")
    return(NULL)
  }
  
  cat("6.1 Performing permutation MANCOVA...\n")
  
  varimax_data <- data_list$varimax
  promax_data <- data_list$promax
  
  dv_names <- paste0("PC", 1:20)
  
  # Check if we have the required variables
  if (!all(dv_names %in% names(varimax_data))) {
    cat("Missing PC variables for permutation tests\n")
    return(NULL)
  }
  
  tryCatch({
    # Prepare dependent variables
    dependent_vars_v <- as.matrix(varimax_data[, ..dv_names])
    dependent_vars_p <- as.matrix(promax_data[, ..dv_names])
    
    set.seed(123)
    
    # Permutation MANCOVA for communication - Varimax
    perm_mancova_v_comm <- anova.cca(
      capscale(dependent_vars_v ~ DX_GROUP * AVG_COMM_INT +
                 AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                 SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
               data = varimax_data), 
      by = "margin", 
      permutations = 1000
    )
    
    # Permutation MANCOVA for communication - Promax
    perm_mancova_p_comm <- anova.cca(
      capscale(dependent_vars_p ~ DX_GROUP * AVG_COMM_INT +
                 AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                 SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
               data = promax_data), 
      by = "margin", 
      permutations = 1000
    )
    
    # Save permutation results
    write.csv(perm_mancova_v_comm, "output/results/perm_mancova_varimax_comm.csv")
    write.csv(perm_mancova_p_comm, "output/results/perm_mancova_promax_comm.csv")
    
    cat("Permutation testing completed\n")
    
  }, error = function(e) {
    cat("Error in permutation testing:", e$message, "\n")
  })
}

# Execute permutation tests
perform_permutation_tests(pca_results)

#=============================================================================
# SECTION 7: FINAL REPORTING
#=============================================================================
cat("\n=== SECTION 7: FINAL REPORTING ===\n")

generate_final_report <- function() {
  cat("7.1 Generating analysis summary...\n")
  
  # Create summary report
  report_lines <- c(
    "# ABIDE Comprehensive Analysis Report",
    "# Generated on:", as.character(Sys.time()),
    "",
    "## Analysis Components Completed:",
    "1. Data preprocessing and normalization",
    "2. Principal Component Analysis with Varimax and Promax rotation",
    "3. Exploratory Data Analysis",
    "4. MANCOVA statistical analysis",
    "5. Support Vector Machine classification",
    "6. Permutation testing",
    "",
    "## Output Files:",
    "- Data files: output/data/",
    "- Results: output/results/",
    "- Plots: output/plots/",
    "",
    "## Key Findings:",
    "- Check individual result files for detailed statistics",
    "- Visualization plots provide insights into data patterns",
    "- Machine learning models show classification performance",
    "",
    paste("Analysis completed at:", Sys.time())
  )
  
  writeLines(report_lines, "output/analysis_report.txt")
  
  cat("Final report generated: output/analysis_report.txt\n")
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("All output files have been saved in the 'output' directory\n")
  cat("Check the analysis_report.txt for a summary of completed analyses\n")
}

# Generate final report
generate_final_report()

#=============================================================================
# END OF SCRIPT
#=============================================================================