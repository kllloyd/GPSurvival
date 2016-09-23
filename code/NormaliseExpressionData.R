NormaliseExpressionData <- function(trainingTestStructure,normaliseFlag,winsoriseFlag){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# To normalised expression data so each gene has mean = 0, sd = 1 																
	# Also can Winsorise the data, i.e. reassign any values outside the (5%,95%) quantiles to these values. This is to remove outliers.
	# Mean, sd and quantiles are calculated on the training data and then applied to both the training and test sets. 			
	#-------------------------------------------------------------------------------------------------------#	

	trainingData 		<- trainingTestStructure$trainingData
	testData 			<- trainingTestStructure$testData
	nTraining 			<- dim(trainingData)[1]
	nTest 				<- dim(testData)[1]
	trainingIndices 	<- 1:nTraining
	testIndices 		<- (nTraining+1):(nTraining+nTest)
	dimension 			<- dim(trainingData)[2]

	trainingMean 		<- apply(trainingData,2,mean)
	trainingSD 			<- apply(trainingData,2,sd)

	dataAll 			<- rbind(trainingData,testData)
	dataAllNormalised 	<- dataAll

	if(normaliseFlag){
		for(g in 1:dimension){
			for(h in 1:(nTraining+nTest)){
				dataAllNormalised[h,g] <- (dataAll[h,g]-trainingMean[g])/trainingSD[g]
				if(is.na(dataAllNormalised[h,g])&trainingSD[g]!=0) cat('dataAll =',dataAll[h,g],', trainingMean =',trainingMean[g],', trainingSD =',trainingSD[g],fill=TRUE)
			}
		}
		dataAllNormalised <- dataAllNormalised[,trainingSD!=0,drop=FALSE]
		dimension <- dim(dataAllNormalised)[2]
	}

	if(winsoriseFlag){
		quantiles 	<- sapply(1:dimension,function(x) quantile(dataAllNormalised[1:nTraining,x],c(0.05,0.95)))
		for(j in 1:dimension){
			for(i in 1:(nTraining+nTest)){
				if(dataAllNormalised[i,j]<quantiles[1,j]) dataAllNormalised[i,j] <- quantiles[1,j]
				if(dataAllNormalised[i,j]>quantiles[2,j]) dataAllNormalised[i,j] <- quantiles[2,j]
			}
		}
	}
	
	trainingTestStructure$trainingData 	<- dataAllNormalised[trainingIndices,,drop=FALSE]
	trainingTestStructure$testData		<- dataAllNormalised[testIndices,,drop=FALSE]
	trainingTestStructure$dimension 	<- dimension

	return(trainingTestStructure)

}