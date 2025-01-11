library(tidyr)
library(readxl)
library(stringr)
library(readr)

###CC400-CPAC-Caltech
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400 - CPAC-Caltech") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
    data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-CMU
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-CMU") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-KKI
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-KKI") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-Leuven1
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-LEUVEN_1") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-Leuven2
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-LEUVEN_2") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-MaxMun
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-MaxMun") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-NYU
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-NYU") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-OHSU
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-OHSU") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-OLIN
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-OLIN") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-OLIN
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-PITT") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-SBL
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-SBL") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-SDSU
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-SDSU") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-STANFORD
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-STANFORD") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-TRINITY
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-TRINITY") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-UCLA1
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-UCLA1") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-UCLA2
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-UCLA2") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-UMI1
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-UMI1") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-USM
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-USM") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

###CC400-CPAC-YALE
setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
path <- setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset]")
folder <- "CC400_CPAC/" ####Change Folder####
X <- read_xlsx("URL_ABIDE_Preprocessing_TS_Data.xlsx", sheet = "CC400-CPAC-YALE") ####Change Sheet####  

list()
URL <- as.data.frame((X$URL))

n <- nrow(URL)
for (i in 1:n) {
  data <- read_tsv(URL[i,1])
  filename <- sub(".*\\/", "", URL[i,1])
  file <-  as.character(paste0(folder,filename,".csv"))
  write.csv(data, file)
}

