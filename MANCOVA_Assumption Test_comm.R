library(MVN)  # For Mardia's test
library(car)  # For Levene's test
library(biotools)  # For Box's M test
library(ggplot2)  # For visualization
library(dplyr)
library(tidyr)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\06.1 MANCOVA Assumption Test]")

# Load datasets
pdata <- read.csv("promax_alldata.csv")
vdata <- read.csv("varimax_alldata.csv")

# Define dependent variables
dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, dv_names]
dependent_vars_v <- vdata[, dv_names]

# Independent variable
pdata$DX_GROUP <- as.factor(pdata$DX_GROUP)
vdata$DX_GROUP <- as.factor(vdata$DX_GROUP)
independent_var <- "DX_GROUP"

# Covariates
covariates <- c("AVG_COMM_INT", "AGE_AT_SCAN", "SEX")

# Variance Inflation Factor (VIF) Test for Multicollinearity
vif_test_p <- vif(lm(PC1 ~ DX_GROUP + AVG_COMM_INT + AGE_AT_SCAN + SEX, data = pdata))
print(vif_test_p)
# DX_GROUP AVG_COMM_INT  AGE_AT_SCAN          SEX 
#1.856257     1.788226     1.057909     1.016888 
vif_test_v <- vif(lm(PC1 ~ DX_GROUP + AVG_COMM_INT + AGE_AT_SCAN + SEX, data = vdata))
print(vif_test_v)
#DX_GROUP AVG_COMM_INT  AGE_AT_SCAN          SEX 
#1.856257     1.788226     1.057909     1.016888 

# Mardia's Test for Multivariate Normality
mardia_result <- mvn(dependent_vars_p, mvnTest = "mardia")
mardia_result_v <- mvn(dependent_vars_v, mvnTest = "mardia")
print(mardia_result$multivariateNormality)
#Test        Statistic p value Result
#1 Mardia Skewness 12416.3208031077       0     NO
#2 Mardia Kurtosis 127.643267471169       0     NO

print(mardia_result_v$multivariateNormality)
#Test        Statistic p value Result
#1 Mardia Skewness  13308.704463159       0     NO
#2 Mardia Kurtosis 138.284671981468       0     NO


# Levene's Test for Homogeneity of Variance
levene_results <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = pdata)
})
names(levene_results) <- dv_names
print(levene_results)

levene_results_v <- lapply(dv_names, function(dv) {
  leveneTest(as.formula(paste(dv, "~", independent_var)), data = vdata)
})
names(levene_results_v) <- dv_names
print(levene_results_v)

# Box's M Test for Homogeneity of Covariance Matrices
box_test <- boxM(dependent_vars_p, pdata[[independent_var]]) #Chi-Sq (approx.) = 646.36, df = 210, p-value < 2.2e-16
box_test_v <- boxM(dependent_vars_v, vdata[[independent_var]]) #Chi-Sq (approx.) = 656.81,df = 210, p-value <2.2e-16
print(box_test)
print(box_test_v)

# Linear Relationship Between Covariates and Dependent Variables
plot_data <- pdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list <- lapply(covariates, function(cov) {
  ggplot(plot_data, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables (Promax)"))
})
print(plot_list)

plot_data_v <- vdata %>%
  pivot_longer(cols = all_of(dv_names), names_to = "Dependent", values_to = "Value")

plot_list_v <- lapply(covariates, function(cov) {
  ggplot(plot_data_v, aes_string(x = cov, y = "Value", color = independent_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~Dependent, scales = "free_y") +
    labs(title = paste("Linear Relationship of", cov, "with Dependent Variables (Varimax)"))
})
print(plot_list_v)

# Testing Equality of Regression Slopes via Interaction Terms
mancova_model <- manova(as.matrix(dependent_vars_p) ~ DX_GROUP * AVG_COMM_INT +
                          AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                          SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                        data = pdata)

mancova_model_v <- manova(as.matrix(dependent_vars_v) ~ DX_GROUP * AVG_COMM_INT +
                            AGE_AT_SCAN * DX_GROUP + AGE_AT_SCAN * AVG_COMM_INT +
                            SEX * DX_GROUP + SEX * AVG_COMM_INT + AGE_AT_SCAN * SEX,
                          data = vdata)

anova_mancova <- summary.aov(mancova_model)
anova_mancova_v <- summary.aov(mancova_model_v)
print(anova_mancova)
print(anova_mancova_v)
