library(dplyr)                                                 
library(data.table)
library(e1071)
library(caret)
library(pROC)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\07 Support Vector Machines]")

# Load the dataset
vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

##varimax dataset training and test 

vdata$DX_GROUP <- as.factor(vdata$DX_GROUP)

# Set seed for reproducibility
set.seed(123)

# Add a unique row identifier to vdata
vdata <- vdata %>%
  mutate(RowID = row_number())

# Create testing dataset with 10% of the data for each level of DX_GROUP 
testing_vdata <- vdata %>%
  filter(DX_GROUP %in% c("1", "2")) %>%
  group_by(DX_GROUP) %>%
  slice_sample(prop = 0.1) %>%
  ungroup()

# Create training dataset by excluding the rows in testing_vdata
training_vdata <- vdata %>%
  filter(!RowID %in% testing_vdata$RowID)

# Drop the RowID column from the final datasets if not needed
testing_vdata <- testing_vdata %>% select(-RowID)
training_vdata <- training_vdata %>% select(-RowID)

# Save files
write.csv(testing_vdata,"testing_vdata.csv")
write.csv(training_vdata, "training_vdata.csv")

##promax dataset training and test 
# Ensure DX_GROUP is a factor
pdata$DX_GROUP <- as.factor(pdata$DX_GROUP)

# Set seed for reproducibility
set.seed(123)

# Add a unique row identifier to pdata
pdata <- pdata %>%
  mutate(RowID = row_number())

# Create testing dataset with 10% of the data for each level of DX_GROUP (1 and 2)
testing_pdata <- pdata %>%
  filter(DX_GROUP %in% c("1", "2")) %>%
  group_by(DX_GROUP) %>%
  slice_sample(prop = 0.1) %>%
  ungroup()

# Create training dataset by excluding the rows in testing_vdata
training_pdata <- pdata %>%
  filter(!RowID %in% testing_pdata$RowID)

# Drop the RowID column from the final datasets if not needed
testing_pdata <- testing_pdata %>% select(-RowID)
training_pdata <- training_pdata %>% select(-RowID)

# Save files
write.csv(testing_pdata,"testing_pdata.csv")
write.csv(training_pdata, "training_pdata.csv")

##SVM implementation
# Ensure DX_GROUP is a factor
training_vdata$DX_GROUP <- as.factor(training_vdata$DX_GROUP)
training_pdata$DX_GROUP <- as.factor(training_pdata$DX_GROUP)

# Rename class levels to valid R variable names
levels(training_vdata$DX_GROUP) <- make.names(levels(training_vdata$DX_GROUP))
levels(training_pdata$DX_GROUP) <- make.names(levels(training_pdata$DX_GROUP))

# Define predictors and dependent variable
predictors <- paste0("PC", 1:20)
dependent_var <- "DX_GROUP"

# Create custom AUC metric for cross-validation
custom_auc <- function(data, lev = NULL, model = NULL) {
  roc_obj <- roc(data$obs, as.numeric(data[, "pred"]))
  auc_value <- auc(roc_obj)
  return(c(AUC = as.numeric(auc_value)))
}

# Define the custom SVM training function
train_svm <- function(training_data) {
  # Ensure DX_GROUP is a factor
  training_data$DX_GROUP <- as.factor(training_data$DX_GROUP)
  
  # Train control with class probabilities enabled
  train_control <- trainControl(
    method = "cv",              
    number = 5,                 
    classProbs = TRUE,          
    savePredictions = TRUE,     # Save predictions for AUC
    summaryFunction = twoClassSummary,  
    verboseIter = TRUE          
  )
  
  # Define grid for hyperparameter tuning
  svm_grid <- expand.grid(
    C = 2^(-5:5),               
    sigma = 2^(-5:5)            
  )
  
  # Train SVM model
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

# Train SVM for training_vdata
cat("Training SVM for training_vdata...\n")
svm_vdata <- train_svm(training_vdata)

# Train SVM for training_pdata
cat("\nTraining SVM for training_pdata...\n")
svm_pdata <- train_svm(training_pdata)

# Display summaries
cat("\nSVM Model Summary for training_vdata:\n")
print(svm_vdata) #sigma = 0.03125 and C = 8

cat("\nSVM Model Summary for training_pdata:\n")
print(svm_pdata) #sigma = 0.03125 and C = 8

# Evaluate the final model AUC for training_vdata
roc_vdata <- roc(
  response = svm_vdata$pred$obs,
  predictor = svm_vdata$pred$X1, 
  levels = rev(levels(svm_vdata$pred$obs)) 
)
cat("\nAUC for training_vdata SVM model:\n")
print(auc(roc_vdata)) #Area under the curve: 0.6915 CV - Insample

# Evaluate the final model AUC for training_pdata
roc_pdata <- roc(
  response = svm_pdata$pred$obs, 
  predictor = svm_pdata$pred$X1, 
  levels = rev(levels(svm_pdata$pred$obs))
)
cat("\nAUC for training_pdata SVM model:\n")
print(auc(roc_pdata)) #Area under the curve: 0.6997 CV - Insample

# Apply the final SVM model to vdata
cat("Applying SVM to training_vdata...\n")
svm_vdata_predictions_all <- predict(svm_vdata, newdata = vdata, type = "prob")
svm_vdata_classes_all <- predict(svm_vdata, newdata = vdata, type = "raw")

# Map predicted labels to match the DX_GROUP levels - vdata
svm_vdata_classes_all <- as.character(svm_vdata_classes_all)  
svm_vdata_classes_all[svm_vdata_classes_all == "X1"] <- "1"
svm_vdata_classes_all[svm_vdata_classes_all == "X2"] <- "2"
svm_vdata_classes_all <- factor(svm_vdata_classes_all, levels = levels(vdata$DX_GROUP))  

# Apply the final SVM model to pdata
cat("Applying SVM to training_pdata...\n")
svm_pdata_predictions_all <- predict(svm_pdata, newdata = pdata, type = "prob")
svm_pdata_classes_all <- predict(svm_pdata, newdata = pdata, type = "raw")

# Map predicted labels to match the DX_GROUP levels - pdata
svm_pdata_classes_all <- as.character(svm_pdata_classes_all)  
svm_pdata_classes_all[svm_pdata_classes_all == "X1"] <- "1"
svm_pdata_classes_all[svm_pdata_classes_all == "X2"] <- "2"
svm_pdata_classes_all <- factor(svm_pdata_classes_all, levels = levels(pdata$DX_GROUP))  

#Combine prediction data to actuals
predictions_vdata_all <- as.data.frame(svm_vdata_classes_all)
predictions_pdata_all <- as.data.frame(svm_pdata_classes_all)
svm_oos_output_vdata_all  <- cbind(vdata, predictions_vdata_all)
svm_oos_output_pdata_all  <- cbind(pdata, predictions_pdata_all)

write.csv(svm_oos_output_vdata_all, "svm_oos_output_vdata_all.csv")
write.csv(svm_oos_output_pdata_all, "svm_oos_output_pdata_all.csv")

##Performance evaluation in OOS data
# Apply the final SVM model to testing_vdata
cat("Applying SVM to testing_vdata...\n")
svm_vdata_predictions <- predict(svm_vdata, newdata = testing_vdata, type = "prob")
svm_vdata_classes <- predict(svm_vdata, newdata = testing_vdata, type = "raw")

# Evaluate the performance on testing_vdata
roc_vdata <- roc(
  testing_vdata$DX_GROUP,
  svm_vdata_predictions[, 2],  
)
auc_vdata <- auc(roc_vdata)
cat("\nAUC for testing_vdata SVM model:\n")
print(auc_vdata) #Area under the curve: 0.8056

# Map predicted labels to match the DX_GROUP levels - vdata
svm_vdata_classes <- as.character(svm_vdata_classes)  
svm_vdata_classes[svm_vdata_classes == "X1"] <- "1"
svm_vdata_classes[svm_vdata_classes == "X2"] <- "2"
svm_vdata_classes <- factor(svm_vdata_classes, levels = levels(testing_vdata$DX_GROUP))  


# Confusion matrix for testing_vdata
testing_vdata$DX_GROUP <- factor(testing_vdata$DX_GROUP, levels = c("1", "2"))
svm_vdata_classes <- factor(svm_vdata_classes, levels = levels(testing_vdata$DX_GROUP))
conf_matrix_vdata <- confusionMatrix(svm_vdata_classes, testing_vdata$DX_GROUP)
cat("\nConfusion Matrix for testing_vdata SVM model:\n")
print(conf_matrix_vdata)

# Apply the final SVM model to testing_pdata
cat("\nApplying SVM to testing_pdata...\n")
svm_pdata_predictions <- predict(svm_pdata, newdata = testing_pdata, type = "prob")
svm_pdata_classes <- predict(svm_pdata, newdata = testing_pdata, type = "raw")

# Map predicted labels to match the DX_GROUP levels - pvdata
svm_pdata_classes <- as.character(svm_pdata_classes)  
svm_pdata_classes[svm_pdata_classes == "X1"] <- "1"
svm_pdata_classes[svm_pdata_classes == "X2"] <- "2"
svm_pdata_classes <- factor(svm_pdata_classes, levels = levels(testing_pdata$DX_GROUP))  

# Evaluate the performance on testing_pdata
roc_pdata <- roc(
testing_pdata$DX_GROUP,
svm_pdata_predictions[, 2],  
levels = rev(levels(testing_pdata$DX_GROUP))
)
auc_pdata <- auc(roc_pdata)
cat("\nAUC for testing_pdata SVM model:\n")
print(auc_pdata) #Area under the curve: 0.7659

# Confusion matrix for testing_pdata
conf_matrix_pdata <- confusionMatrix(svm_pdata_classes, testing_pdata$DX_GROUP)
cat("\nConfusion Matrix for testing_pdata SVM model:\n")
print(conf_matrix_pdata)


# Extract confusion matrix details into data frames

# For conf_matrix_pdata
conf_matrix_pdata_df <- as.data.frame(as.table(conf_matrix_pdata$table))  # Extract table
conf_matrix_pdata_metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Value = c(conf_matrix_pdata$overall["Accuracy"], 
            conf_matrix_pdata$byClass["Sensitivity"], 
            conf_matrix_pdata$byClass["Specificity"])
)

# Combine table and metrics for pdata
conf_matrix_pdata_combined <- list(
  Confusion_Table = conf_matrix_pdata_df,
  Metrics = conf_matrix_pdata_metrics
)

# For conf_matrix_vdata
conf_matrix_vdata_df <- as.data.frame(as.table(conf_matrix_vdata$table))  # Extract table
conf_matrix_vdata_metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Value = c(conf_matrix_vdata$overall["Accuracy"], 
            conf_matrix_vdata$byClass["Sensitivity"], 
            conf_matrix_vdata$byClass["Specificity"])
)

# Combine table and metrics for vdata
conf_matrix_vdata_combined <- list(
  Confusion_Table = conf_matrix_vdata_df,
  Metrics = conf_matrix_vdata_metrics
)

# Save to CSV if needed
write.csv(conf_matrix_pdata_df, "conf_matrix_pdata_table.csv", row.names = FALSE)
write.csv(conf_matrix_pdata_metrics, "conf_matrix_pdata_metrics.csv", row.names = FALSE)

write.csv(conf_matrix_vdata_df, "conf_matrix_vdata_table.csv", row.names = FALSE)
write.csv(conf_matrix_vdata_metrics, "conf_matrix_vdata_metrics.csv", row.names = FALSE)

#Combine prediction data to actuals
predictions_vdata <- as.data.frame(svm_vdata_classes)
predictions_pdata <- as.data.frame(svm_pdata_classes)
svm_oos_output_vdata <- cbind(testing_vdata, predictions_vdata)
svm_oos_output_pdata <- cbind(testing_pdata, predictions_pdata)

write.csv(svm_oos_output_vdata, "svm_oos_output_vdata.csv")
write.csv(svm_oos_output_pdata, "svm_oos_output_pdata.csv")


###IN-SAMPLE ROC CURVES (from cross-validation predictions)

# Set plotting layout: 2 rows, 2 columns
par(mfrow = c(2, 2))  # 2x2 layout for four plots

### IN-SAMPLE ROC CURVES ====

# Varimax - In-sample ROC
roc_cv_vdata <- roc(
  response = svm_vdata$pred$obs,
  predictor = svm_vdata$pred$X1,
  levels = rev(levels(svm_vdata$pred$obs))
)
plot(roc_cv_vdata, col = "blue", main = "In-Sample ROC - Varimax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_cv_vdata), 4)), col = "blue", lwd = 2)

# Promax - In-sample ROC
roc_cv_pdata <- roc(
  response = svm_pdata$pred$obs,
  predictor = svm_pdata$pred$X1,
  levels = rev(levels(svm_pdata$pred$obs))
)
plot(roc_cv_pdata, col = "darkgreen", main = "In-Sample ROC - Promax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_cv_pdata), 4)), col = "darkgreen", lwd = 2)

###OUT-OF-SAMPLE ROC CURVES

# Varimax - Out-of-sample ROC
roc_oos_vdata <- roc(
  response = testing_vdata$DX_GROUP,
  predictor = svm_vdata_predictions[, 2],
  levels = rev(levels(testing_vdata$DX_GROUP))
)
plot(roc_oos_vdata, col = "red", main = "Out-of-Sample ROC - Varimax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_oos_vdata), 4)), col = "red", lwd = 2)

# Promax - Out-of-sample ROC
roc_oos_pdata <- roc(
  response = testing_pdata$DX_GROUP,
  predictor = svm_pdata_predictions[, 2],
  levels = rev(levels(testing_pdata$DX_GROUP))
)
plot(roc_oos_pdata, col = "purple", main = "Out-of-Sample ROC - Promax", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_oos_pdata), 4)), col = "purple", lwd = 2)
