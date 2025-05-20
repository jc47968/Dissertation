library(car)
library(data.table)

vdata <- fread("varimax_alldata.csv")
pdata <- fread("promax_alldata.csv")

# Select only the predictor variables (PC1 to PC20)
predictor_vars <- paste0("PC", 1:20)

# Ensure all predictor variables are numeric
vdata[, (predictor_vars) := lapply(.SD, as.numeric), .SDcols = predictor_vars]
pdata[, (predictor_vars) := lapply(.SD, as.numeric), .SDcols = predictor_vars]

# Fit a linear model with DX_GROUP as the dependent variable and PC1 to PC20 as predictors
lm_vdata <- lm(DX_GROUP ~ ., data = vdata[, c("DX_GROUP", ..predictor_vars), with = FALSE])
lm_pdata <- lm(DX_GROUP ~ ., data = pdata[, c("DX_GROUP", ..predictor_vars), with = FALSE])

# Compute VIF values
vif_vdata <- vif(lm_vdata)
vif_pdata <- vif(lm_pdata)

# Print results
cat("VIF Results for Varimax Dataset:\n")
print(vif_vdata)

cat("\nVIF Results for Promax Dataset:\n")
print(vif_pdata)

# Save VIF results to CSV
write.csv(data.frame(Variable = names(vif_vdata), VIF = vif_vdata), "vif_results_varimax.csv", row.names = FALSE)
write.csv(data.frame(Variable = names(vif_pdata), VIF = vif_pdata), "vif_results_promax.csv", row.names = FALSE)
