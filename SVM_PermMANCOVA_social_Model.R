# Permutation MANCOVA for Social Cognition
library(dplyr)                                                 
library(data.table)
library(car)
library(vegan)  # For permutation tests

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\08 SVM MANCOVA]")

vdata <- fread("svm_oos_output_vdata_all.csv")
pdata <- fread("svm_oos_output_pdata_all.csv")

colnames(vdata)[colnames(vdata) == "svm_vdata_classes_all"] <- "pred"
colnames(pdata)[colnames(pdata) == "svm_pdata_classes_all"] <- "pred"

pdata$pred <- as.factor(pdata$pred)
vdata$pred <- as.factor(vdata$pred)

# Check the number of levels in predictor
if (length(levels(pdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in pdata.")
if (length(levels(vdata$pred)) < 2) stop("Error: Predictor variable 'pred' has less than two levels in vdata.")

# Define the dependent variables (PC1 to PC20)
dv_names <- paste0("PC", 1:20)

dependent_vars_p <- pdata[, ..dv_names]
dependent_vars_v <- vdata[, ..dv_names]

set.seed(123)  

# Permutation MANCOVA for pdata
perm_mancova_p <- anova.cca(capscale(as.matrix(dependent_vars_p) ~ pred * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * pred + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = pdata), by = "margin", permutations = 10000)
print(perm_mancova_p)
write.csv(perm_mancova_p, "PermMANCOVA_svmdata_promax_social.csv")

# Permutation MANCOVA for vdata
perm_mancova_v <- anova.cca(capscale(as.matrix(dependent_vars_v) ~ pred * AVG_SOCIAL_INT +
                                       AGE_AT_SCAN * pred + AGE_AT_SCAN * AVG_SOCIAL_INT +
                                       SEX * pred + SEX * AVG_SOCIAL_INT + AGE_AT_SCAN * SEX,
                                     data = vdata), by = "margin", permutations = 10000)
print(perm_mancova_v)
write.csv(perm_mancova_v, "PermMANCOVA_svmdata_varimax_social.csv")

# Bonferroni Post Hoc Test for Interaction pred:AVG_SOCIAL_INT
bonferroni_interaction_tests <- function(data, predictor, interaction_term) {
  results <- list()
  
  for (dv in dv_names) {
    formula <- as.formula(paste(dv, "~", predictor, "*", interaction_term))
    bonferroni_result <- tryCatch({
      pairwise.t.test(data[[dv]], interaction(data[[predictor]], data[[interaction_term]]), p.adjust.method = "bonferroni")
    }, error = function(e) {
      message(paste("Bonferroni test failed for", dv, "due to insufficient group levels"))
      return(NULL)
    })
    
    if (!is.null(bonferroni_result)) {
      results[[paste(dv, "Bonferroni Interaction")]] <- bonferroni_result$p.value
    }
  }
  return(results)
}

# Perform Bonferroni post hoc tests for pred:AVG_SOCIAL_INT
posthoc_p_bonferroni_interaction <- bonferroni_interaction_tests(pdata, "pred", "AVG_SOCIAL_INT")
posthoc_v_bonferroni_interaction <- bonferroni_interaction_tests(vdata, "pred", "AVG_SOCIAL_INT")

# Print results
print(posthoc_p_bonferroni_interaction)
print(posthoc_v_bonferroni_interaction)
