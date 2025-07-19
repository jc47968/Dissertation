library(dplyr)
library(data.table)
library(multcomp)
library(tidyr)
library(rstatix)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\08 SVM MANCOVA]")

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


# Function to conduct Games-Howell test
run_games_howell <- function(data, factor_var, dataset_name) {
  results <- list()
  
  for (dv in pc_columns) {
    
    # Check if columns exist before proceeding
    required_cols <- c(dv, "pred", factor_var)
    if (!all(required_cols %in% names(data))) {
      message(paste("Skipping", dv, "- missing required columns."))
      next
    }
    
    # Subset and drop NAs
    subset_data <- data[, ..required_cols]
    subset_data <- subset_data[complete.cases(subset_data), ]
    
    # Create interaction group
    subset_data$group <- interaction(subset_data$pred, subset_data[[factor_var]], drop = TRUE)
    
    # Filter to keep only groups with at least 2 observations
    group_counts <- table(subset_data$group)
    valid_groups <- names(group_counts[group_counts >= 2])
    subset_data <- subset_data[subset_data$group %in% valid_groups, ]
    
    # Proceed only if we have at least 2 valid groups
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

# Run the test
gh_pdata_social <- run_games_howell(pdata, "AVG_SOCIAL_INT", "pdata")
gh_vdata_social <- run_games_howell(vdata, "AVG_SOCIAL_INT", "vdata")

gh_pdata_comm <- run_games_howell(pdata, "AVG_COMM_INT", "pdata")
gh_vdata_comm <- run_games_howell(vdata, "AVG_COMM_INT", "vdata")

# Combine and save
final_social <- bind_rows(gh_pdata_social, gh_vdata_social)
final_comm <- bind_rows(gh_pdata_comm, gh_vdata_comm)

write.csv(final_social, "Games_Howell_Results_Social.csv", row.names = FALSE)
write.csv(final_comm, "Games_Howell_Results_Comm.csv", row.names = FALSE)
