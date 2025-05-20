library(dplyr)
library(data.table)
library(multcomp)
library(tidyr)

pdata <- fread("svm_oos_output_pdata_all.csv")
vdata <- fread("svm_oos_output_vdata_all.csv")

# Select relevant columns (PC1 to PC20)
pc_columns <- paste0("PC", 1:20)

pdata$pred <- as.factor(pdata$svm_pdata_classes_all)
vdata$pred <- as.factor(vdata$svm_vdata_classes_all)


#Tukey's HSD for AVG_SOCIAL_INT
pdata$AVG_SOCIAL_INT <- as.factor(pdata$AVG_SOCIAL_INT)
vdata$AVG_SOCIAL_INT <- as.factor(vdata$AVG_SOCIAL_INT)

# Function to conduct Tukey's HSD test
run_tukey_hsd <- function(data, dataset_name) {
  results <- list()
  
  for (dv in pc_columns) {
    if (!dv %in% colnames(data)) {
      message(paste("Skipping", dv, "because it is not found in the dataset."))
      next
    }
    
    formula <- as.formula(paste(dv, "~ pred * AVG_SOCIAL_INT"))
    
    # Fit ANOVA model
    model <- tryCatch({
      aov(formula, data = data)
    }, error = function(e) {
      message(paste("ANOVA failed for", dv, ":", e$message))
      return(NULL)
    })
    
    # Perform Tukey's HSD test
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
  
  # Convert results to tidy format
  tidy_results <- list()
  
  for (dv in names(results)) {
    if (!is.null(results[[dv]])) {
      tidy_data <- as.data.frame(results[[dv]][[1]]) %>%
        tibble::rownames_to_column(var = "Comparison") %>%
        mutate(Dependent_Variable = dv)
      tidy_results[[dv]] <- tidy_data
    }
  }
  
  # Combine results into a single data frame
  final_results <- bind_rows(tidy_results)
  final_results$Dataset <- dataset_name
  
  return(final_results)
}

# Run Tukey's HSD for both datasets
tukey_results_pdata <- run_tukey_hsd(pdata, "pdata")
tukey_results_vdata <- run_tukey_hsd(vdata, "vdata")

# Combine results
combined_results <- bind_rows(tukey_results_pdata, tukey_results_vdata)

# Display and save results
print(combined_results)
write.csv(combined_results, "Tukey_HSD_Results_Social.csv", row.names = FALSE)

#Tukey's HSD for AVG_COMM_INT
pdata$AVG_COMM_INT <- as.factor(pdata$AVG_COMM_INT)
vdata$AVG_COMM_INT <- as.factor(vdata$AVG_COMM_INT)

# Function to conduct Tukey's HSD test
run_tukey_hsd <- function(data, dataset_name) {
  results <- list()
  
  for (dv in pc_columns) {
    if (!dv %in% colnames(data)) {
      message(paste("Skipping", dv, "because it is not found in the dataset."))
      next
    }
    
    formula <- as.formula(paste(dv, "~ pred * AVG_COMM_INT"))
    
    # Fit ANOVA model
    model <- tryCatch({
      aov(formula, data = data)
    }, error = function(e) {
      message(paste("ANOVA failed for", dv, ":", e$message))
      return(NULL)
    })
    
    # Perform Tukey's HSD test
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
  
  # Convert results to tidy format
  tidy_results <- list()
  
  for (dv in names(results)) {
    if (!is.null(results[[dv]])) {
      tidy_data <- as.data.frame(results[[dv]][[1]]) %>%
        tibble::rownames_to_column(var = "Comparison") %>%
        mutate(Dependent_Variable = dv)
      tidy_results[[dv]] <- tidy_data
    }
  }
  
  # Combine results into a single data frame
  final_results <- bind_rows(tidy_results)
  final_results$Dataset <- dataset_name
  
  return(final_results)
}

# Run Tukey's HSD for both datasets
tukey_results_pdata <- run_tukey_hsd(pdata, "pdata")
tukey_results_vdata <- run_tukey_hsd(vdata, "vdata")

# Combine results
combined_results <- bind_rows(tukey_results_pdata, tukey_results_vdata)

# Display and save results
print(combined_results)
write.csv(combined_results, "Tukey_HSD_Results_Comm.csv", row.names = FALSE)
