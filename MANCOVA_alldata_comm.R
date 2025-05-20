library(dplyr)                                                 
library(data.table)
library(car)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\06 MANCOVA_Alldata]")

# Load the dataset
vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

# Define the dependent variables (PC1 to PC20)
dependent_vars <- paste0("PC", 1:20)

# Create the formula for MANCOVA
mancova_formula <- as.formula(
  paste("cbind(", paste(dependent_vars, collapse = ", "), 
        ") ~ DX_GROUP * AVG_COMM_INT + AGE_AT_SCAN*DX_GROUP + AGE_AT_SCAN*AVG_COMM_INT +
        SEX*DX_GROUP + SEX*AVG_COMM_INT + AGE_AT_SCAN * SEX")
)

# Perform MANCOVA for varimax
varimax_mancova_result <- manova(mancova_formula, data = vdata)

# View the MANCOVA summary for varimax
varimax_summary_mancova <- summary(varimax_mancova_result, test = "Wilks")

# Detailed summary for each dependent variable for varimax
varimax_summary_aov <- summary.aov(varimax_mancova_result)

# Perform MANCOVA for promax
promax_mancova_result <- manova(mancova_formula, data = pdata)

# View the MANCOVA summary for promax
promax_summary_mancova <- summary(promax_mancova_result, test = "Wilks")

# Detailed summary for each dependent variable for promax 
promax_summary_aov <- summary.aov(promax_mancova_result)
print(promax_summary_aov)

# Convert MANCOVA summary to a data frame - varimax
varimax_summary_df <- as.data.frame(varimax_summary_mancova$stats)
varimax_summary_df <- cbind(Variable = rownames(varimax_summary_df), varimax_summary_df)
rownames(varimax_summary_df) <- NULL

# Convert MANCOVA summary to a data frame - promax
promax_summary_df <- as.data.frame(promax_summary_mancova$stats)
promax_summary_df <- cbind(Variable = rownames(promax_summary_df), promax_summary_df)
rownames(promax_summary_df) <- NULL


# Convert ANOVA summary to a data frame - varimax
varimax_anova_list <- lapply(varimax_summary_aov, function(x) {
  as.data.frame(x)
})

# Convert ANOVA summary to a data frame - promax
promax_anova_list <- lapply(promax_summary_aov, function(x) {
  as.data.frame(x)
})

# Combine into a single data frame - varimax
varimax_detailed_summary_df <- do.call(rbind, lapply(names(varimax_anova_list), function(var) {
  df <- varimax_anova_list[[var]]
  df$Dependent_Variable <- var  
  return(df)
}))

# Combine into a single data frame - promax
promax_detailed_summary_df <- do.call(rbind, lapply(names(promax_anova_list), function(var) {
  df <- promax_anova_list[[var]]
  df$Dependent_Variable <- var  
  return(df)
}))

#Save files
write.csv(varimax_detailed_summary_df, "varimax_detailed_summary_df[comm].csv")
write.csv(varimax_summary_df, "varimax_summary_df[comm].csv")
write.csv(promax_detailed_summary_df, "promax_detailed_summary_df[comm].csv")
write.csv(promax_summary_df, "promax_summary_df[comm].csv")


