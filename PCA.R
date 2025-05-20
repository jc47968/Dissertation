library(GPArotation)
library(psych)
library(dplyr)                                                 
library(data.table)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\04 Princple Component Analysis]")

merged_data <- fread("merged_data.csv")

# Start the runtime meter
start_time <- Sys.time()

# Select the variables for PCA 
# Get the column range between "X.10_X.1" and "X.99_X.98"
start_col <- which(names(merged_data) == "X.10_X.1")
end_col <- which(names(merged_data) == "X.99_X.98")

# Subset columns in the range
selected_vars <- merged_data[,start_col:end_col]
numeric_vars <- na.omit(as.data.frame(lapply(selected_vars, as.numeric)))

# Perform PCA using prcomp
pca_result <- prcomp(numeric_vars, center = TRUE, scale. = TRUE)

# Eigenvalues (variances explained by each component)
eigenvalues <- pca_result$sdev^2

# Save the scree plot
png("scree_plot_prcomp_25_intervals.png", width = 800, height = 600)
num_components <- length(eigenvalues)  # Total number of components
x_intervals <- seq(1, num_components, length.out = 25)  # Create 25 evenly spaced intervals

plot(eigenvalues, type = "b", main = "Scree Plot", xlab = "Principal Component", ylab = "Eigenvalue",
     pch = 19, col = "blue", xaxt = "n")  # Suppress x-axis temporarily
axis(1, at = x_intervals
     , labels = round(x_intervals))  # Add custom x-axis with 25 intervals
abline(h = 1, col = "red", lty = 2)  # Reference line at eigenvalue = 1
dev.off()
cat("Scree plot.png\n")

# Perform PCA loadings
loadings <- as.matrix(pca_result$rotation)

# Extract the loadings for the first 50 components based on scree plot
loadings_50 <- as.matrix(pca_result$rotation[, 1:50])

# Perform Varimax rotation
varimax_result <- varimax(loadings_50)  
rotated_loadings <- varimax_result$loadings

# Function to find the most important features for each component
get_top_features <- function(rotated_loadings, top_n = 5) {
  rotated_loadings <- as.matrix(rotated_loadings)
  if (ncol(rotated_loadings) == 0 || nrow(rotated_loadings) == 0) {
    stop("Error: rotated_loadings has zero rows or columns.")
  }
    # Initialize a list to store results
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

# Get the most important features for Varimax rotated loadings
rotated_loadings_matrix <- as.matrix(rotated_loadings)  # Convert to matrix
top_features_varimax <- get_top_features(rotated_loadings_matrix, top_n = 10)

# Directory to save the CSV files
output_directory <- r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\04 Princple Component Analysis\Varimax Factor Importance]"  
if (!dir.exists(output_directory)) {
  dir.create(output_directory)  # Create the directory if it doesn't exist
}

# Loop through each component and save its loadings to a CSV file
for (i in seq_len(ncol(rotated_loadings))) {
  component_loadings <- rotated_loadings[, i, drop = FALSE]
  filename <- paste0(output_directory, "/component_", i, "_loadings.csv")
  write.csv(component_loadings, filename, row.names = TRUE)
  cat("Saved loadings for Component", i, "to", filename, "\n")
}

# Capture all the rotated loadings into a single data frame
rotated_loadings <- as.matrix(rotated_loadings)

# Restrict to the first 20 components
num_components <- min(20, ncol(rotated_loadings))
selected_loadings <- rotated_loadings[, 1:num_components]

# Convert to a data frame
loadings_dataframe <- as.data.frame(selected_loadings)

# Add feature names as a column for better readability
loadings_dataframe$Feature <- rownames(rotated_loadings)

# Reorder the columns so 'Feature' comes first
loadings_dataframe <- loadings_dataframe[, c("Feature", colnames(loadings_dataframe)[-ncol(loadings_dataframe)])]

# Save the data frame to a CSV file
write.csv(loadings_dataframe, "rotated_loadings_20_components_varimax.csv", row.names = FALSE)

# Create a data frame to capture the highest and lowest values for each variable
extreme_values <- data.frame(
  Feature = rownames(selected_loadings),
  Highest_Component = NA,
  Highest_Value = NA,
  Lowest_Component = NA,
  Lowest_Value = NA
)

# Loop through each variable and find the highest and lowest values
for (i in 1:nrow(selected_loadings)) {
  variable_loadings <- selected_loadings[i, ]
  extreme_values$Highest_Component[i] <- which.max(variable_loadings)
  extreme_values$Highest_Value[i] <- max(variable_loadings)
  extreme_values$Lowest_Component[i] <- which.min(variable_loadings)
  extreme_values$Lowest_Value[i] <- min(variable_loadings)
}

# Save the results to a CSV file
write.csv(extreme_values, "Factor Importance_varimax.csv", row.names = FALSE)

# Calculate varimax rotated scores
numeric_vars_matrix <- as.matrix(numeric_vars)  # Ensure data is numeric
varimax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings)

# Perform Promax rotation
promax_result <- promax(loadings_50)  # Apply Promax rotation
rotated_loadings <- promax_result$loadings

# Function to find the most important features for each component
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
    
    # Sort by absolute loading value
    component_features <- component_features[order(abs(component_features$Loading), decreasing = TRUE), ]
    
    # Handle case where top_n exceeds number of features
    top_n <- min(top_n, nrow(component_features))
    
    # Select the top N features for the component
    top_features <- component_features[1:top_n, ]
    top_features_list[[paste0("Component_", i)]] <- top_features
  }
  
  return(top_features_list)
}
# Get the most important features for Promax rotated loadings
rotated_loadings_matrix <- as.matrix(rotated_loadings)  
top_features_promax<- get_top_features(rotated_loadings_matrix, top_n = 10)

# Directory to save the CSV files
output_directory <- r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\04 Princple Component Analysis\Promax Factor Importance]"  
if (!dir.exists(output_directory)) {
  dir.create(output_directory)  # Create the directory if it doesn't exist
}

# Loop through each component and save its loadings to a CSV file
for (i in seq_len(ncol(rotated_loadings))) {
  component_loadings <- rotated_loadings[, i, drop = FALSE]
    # Create a filename for the component
  filename <- paste0(output_directory, "/component_", i, "_loadings.csv")
    # Save the loadings to a CSV file
  write.csv(component_loadings, filename, row.names = TRUE)
    cat("Saved loadings for Component", i, "to", filename, "\n")
}

# Capture all the rotated loadings into a single data frame
# Ensure `rotated_loadings` is a matrix
rotated_loadings <- as.matrix(rotated_loadings)

# Restrict to the first 20 components
num_components <- min(20, ncol(rotated_loadings))
selected_loadings <- rotated_loadings[, 1:num_components]

# Convert to a data frame
loadings_dataframe <- as.data.frame(selected_loadings)

# Add feature names as a column 
loadings_dataframe$Feature <- rownames(rotated_loadings)

# Reorder the columns so 'Feature' comes first
loadings_dataframe <- loadings_dataframe[, c("Feature", colnames(loadings_dataframe)[-ncol(loadings_dataframe)])]

# Save the data frame to a CSV file
write.csv(loadings_dataframe, "rotated_loadings_20_components_promax.csv", row.names = FALSE)

# Create a data frame to capture the highest and lowest values for each variable
extreme_values <- data.frame(
  Feature = rownames(selected_loadings),
  Highest_Component = NA,
  Highest_Value = NA,
  Lowest_Component = NA,
  Lowest_Value = NA
)

# Loop through each variable and find the highest and lowest values
for (i in 1:nrow(selected_loadings)) {
  # Get the loadings for the current variable across all components
  variable_loadings <- selected_loadings[i, ]
  
  # Find the highest and lowest values and their corresponding components
  extreme_values$Highest_Component[i] <- which.max(variable_loadings)
  extreme_values$Highest_Value[i] <- max(variable_loadings)
  extreme_values$Lowest_Component[i] <- which.min(variable_loadings)
  extreme_values$Lowest_Value[i] <- min(variable_loadings)
}

# Save the results to a CSV file
write.csv(extreme_values, "Factor Importance_promax.csv", row.names = FALSE)

# Calculate promax rotated scores
numeric_vars_matrix <- as.matrix(numeric_vars)  
promax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings)

# Extract rotated loadings
rotated_loadings_varimax <- varimax_result$loadings
rotated_loadings_promax <- promax_result$loadings

# Define a threshold for loading significance
threshold <- 0.00  

# Function to map variables to components
get_variables_by_component_wide <- function(rotated_loadings, threshold) {
  results <- list()
  for (i in 1:ncol(rotated_loadings)) {
    component_vars <- rownames(rotated_loadings)[which(abs(rotated_loadings[, i]) > threshold)]
    results[[paste0("Component_", i)]] <- component_vars
  }
  # Create wide-format data frame
  max_length <- max(sapply(results, length))  
  wide_data <- as.data.frame(do.call(cbind, lapply(results, function(x) {
    c(x, rep(NA, max_length - length(x)))  
  })))
  colnames(wide_data) <- names(results)
  return(wide_data)
}

# Get wide-format data frames for Varimax and Promax
wide_varimax <- get_variables_by_component_wide(rotated_loadings_varimax, threshold)
wide_promax <- get_variables_by_component_wide(rotated_loadings_promax, threshold)

# View the results
cat("\nVarimax Wide-Format Results:\n")
print(wide_varimax)

cat("\nPromax Wide-Format Results:\n")
print(wide_promax)

#Save outputs
write.csv(varimax_scores, "varimax_scores.csv")
write.csv(promax_scores, "promax_scores.csv")
write.csv(wide_varimax, "varimax_pc_variables.csv")
write.csv(wide_promax, "promax_pc_variables.csv")

# Total runtime
total_end_time <- Sys.time()
cat("Total runtime:", round(difftime(total_end_time, start_time, units = "secs"), "seconds\n"))

