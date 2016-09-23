ImputeMissingData <- function(allData,imputeType){
	#---------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #---------------------------------------------------------------------------#
	nRow <- dim(allData)[1]
	nCol <- dim(allData)[2]

	if(imputeType=='ImputeMean'){
		meanValue <- sapply(1:nCol,function(x) mean(allData[,x],na.rm=TRUE))
		for(i in 1:nCol){
			allData[is.na(allData[,i]),i] <- meanValue[i]
		}
	} else if(imputeType=='ImputeMICE'){
		test <- list()
		test[[1]] <- allData
		i <- 1
		while(any(is.na(test[[i]]))){
			imp <- mice(test[[i]],m=5,printFlag=FALSE)
			test[[i+1]] <- complete(imp,1)
			i <- i+1
		}
		allData <- test[[length(test)]]
	} else if (imputeType=='ImputeBioc'){
		allData <- as.data.frame(impute.knn(as.matrix(allData))$data)
	}

	allFeatureNames <- names(allData)

	toReturn <- list('allData'=allData,'allFeatureNames'=allFeatureNames)
	return(toReturn)
}