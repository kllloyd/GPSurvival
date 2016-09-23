# Yuan Yuan 2012-12-10
# This script read the tab-separated file and handles proper rownames/colnames assignment, also deals with the end-line tab
# also remove the non-informative features
## EDITED K Lloyd 2016_04_29 ##
## loadEntity is depreciated ##
library(randomForest)

synapse.read <- function(synapse.ID)
  {
    synapse.data <- synGet(synapse.ID)
    fileName <- synapse.data@filePath
    return(fileName)
  }

myRead <- function(synapse.ID)
  {
    fileName <- synapse.read(synapse.ID)
    
    data <- read.table(fileName, header=T,na.strings=c("[Pending]","[Not Available]","[Not Applicable]","null","null ","NA"),  sep="\t", quote="")
    samples <- as.character(data[,1])
    data <- data[,2:(ncol(data)-1)]
    rownames(data) <- samples
    print(paste(fileName, ": ", nrow(data), "(samples),", ncol(data), "(features)"))
    data <- na.roughfix(data) # impute the missing value in x, require randomForest library
    valid.cols <- which(apply(data, 2, sd)> 0) # remove the non-informative ones, with same values (e.g., 0) across all samples
    print(paste((ncol(data)-length(valid.cols)), "invalid cols were removed."))
    data <- data[,valid.cols]    
    return (data)
  }

# deals with clinical data, in which factors and numbers may be mixed but not check for flat values
myRead.simple <- function(synapse.ID)
  {
    fileName <- synapse.read(synapse.ID)
    
    data <- read.table(fileName, header=T, na.strings=c("[Pending]","[Not Available]","[Not Applicable]","null","null ","NA"), sep="\t", quote="")
    samples <- as.character(data[,1])
    data <- data[,2:(ncol(data)-1)]
    rownames(data) <- samples
    print(paste(fileName, ": ", nrow(data), "(samples),", ncol(data), "(features)"))
    data <- na.roughfix(data) # impute the missing value in x, require randomForest library    
    return (data)
  }
