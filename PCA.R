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
loadings <- as.matrix(pca_result$rotation)  # Extract loadings

# Extract the loadings for the first 50 components based on scree plot
loadings_50 <- as.matrix(pca_result$rotation[, 1:50])

# Perform Varimax rotation
varimax_result <- varimax(loadings_50)  # Apply Varimax rotation
rotated_loadings <- varimax_result$loadings
numeric_vars_matrix <- as.matrix(numeric_vars)  # Ensure data is numeric
varimax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings)

# Perform Promax rotation
promax_result <- promax(loadings_50)  # Apply Promax rotation
rotated_loadings <- promax_result$loadings
# Calculate Varimax rotated scores
numeric_vars_matrix <- as.matrix(numeric_vars)  # Ensure data is numeric
promax_scores <- as.data.frame(numeric_vars_matrix %*% rotated_loadings)

# Extract rotated loadings
rotated_loadings_varimax <- varimax_result$loadings
rotated_loadings_promax <- promax_result$loadings

# Define a threshold for loading significance
threshold <- 0.00  # Set threshold

# Function to map variables to components
get_variables_by_component_wide <- function(rotated_loadings, threshold) {
  results <- list()
  for (i in 1:ncol(rotated_loadings)) {
    component_vars <- rownames(rotated_loadings)[which(abs(rotated_loadings[, i]) > threshold)]
    results[[paste0("Component_", i)]] <- component_vars
  }
  # Create wide-format data frame
  max_length <- max(sapply(results, length))  # Find the maximum number of variables in any component
  wide_data <- as.data.frame(do.call(cbind, lapply(results, function(x) {
    c(x, rep(NA, max_length - length(x)))  # Fill shorter columns with NA
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

