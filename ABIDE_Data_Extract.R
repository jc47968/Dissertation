library(readxl)
library(readr)
library(stringr)

# Set working directory
setwd("C:/Users/jason/OneDrive/Desktop/NCU/Dissertation Dataset")
folder <- "CC400_CPAC/"

# Get all sheet names
sheets <- excel_sheets("URL_ABIDE_Preprocessing_TS_Data.xlsx")

# Loop through each sheet
for (sheet in sheets) {
  cat("Processing:", sheet, "\n")
  
  # Read the sheet
  df <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = sheet)
  
  # Extract URL column
  URLs <- as.data.frame(df$URL)
  
  # Download and save each file
  for (i in 1:nrow(URLs)) {
    url <- URLs[i, 1]
    filename <- sub(".*\\/", "", url)
    filepath <- paste0(folder, filename, ".csv")
    
    data <- read_tsv(url)
    write.csv(data, filepath)
  }
}
