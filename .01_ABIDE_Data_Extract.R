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
