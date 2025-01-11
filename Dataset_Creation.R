library(dplyr)                                                 
library(data.table)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\05 Dataset Creation]")
merged_data <- fread("merged_data.csv")
varimax_scores <- fread("varimax_scores.csv")
promax_scores <- fread("promax_scores.csv")

#Create datasets
alldata <- merged_data %>% select(SUB_ID, SITE_ID, DX_GROUP, DSM_IV_TR,
                                  AGE_AT_SCAN, SEX, AVG_SOCIAL, AVG_COMM)

# Convert AVG_SOCIAL and AVG_COMM to integers between 0 and 10
alldata$AVG_SOCIAL_INT <- round(alldata$AVG_SOCIAL * 10)
alldata$AVG_COMM_INT <- round(alldata$AVG_COMM * 10)

#Combine datset to PCA data
varimax_alldata <- combined_data <- cbind(alldata, varimax_scores)
varimax_alldata <- varimax_alldata %>% select (-V1)
promax_alldata <- combined_data <- cbind(alldata, promax_scores)
promax_alldata <- promax_alldata %>% select (-V1)

#Save data
write.csv(varimax_alldata, "varimax_alldata.csv")
write.csv(promax_alldata, "promax_alldata.csv")


