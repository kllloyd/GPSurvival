GenerateData <- function(dataOptionsStructure,outerFolder,nReps){
	#-----------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#-----------------------------------------------------------------------------------------------------#
	# Function applies MakeSyntheticData.R to generate data relevant to the problem.
	# If data requirements change, change this function (e.g. if test set needs to be random)
	# If data parameters change, change them in runGPLaplaceApproximation.R
	#-----------------------------------------------------------------------------------------------------#

	dataSource 		<- dataOptionsStructure$dataSource
	censoringType 	<- dataOptionsStructure$censoringType
	nTraining 		<- dataOptionsStructure$nTraining
	nTest 			<- dataOptionsStructure$nTest
	dimension 		<- dataOptionsStructure$dimension

	if(dataSource=='Generate'){
		dataStructure 		<- list()
		allX 				<- list()
		allY 				<- list()
		allYPreCensoring 	<- list()
		censoredStructure 	<- list()
		allEvents 			<- list()
		for(i in 1:nReps){
			dataOptionsStructure$nSamples 		<- nTraining + nTest
			dataStructure[[i]] 					<- MakeSyntheticData(dataOptionsStructure)
			allX[[i]] 							<- dataStructure[[i]]$data
			allY[[i]] 							<- dataStructure[[i]]$targets
			if(censoringType!='None'){
				allYPreCensoring[[i]] 			<- dataStructure[[i]]$targets
				censoredStructure[[i]] 			<- CensorData(dataStructure[[i]],dataOptionsStructure)
				allY[[i]] 						<- censoredStructure[[i]]$targets
				allEvents[[i]] 					<- censoredStructure[[i]]$events
			}
		}
	} else if(dataSource=='CuratedOvarian'){
		dataSet 					<- dataOptionsStructure$dataSet
		geneExpressionFlag 			<- dataOptionsStructure$geneExpressionFlag
		clinicalFeaturesFlag 		<- dataOptionsStructure$clinicalFeaturesFlag
		chosenClinicalFeatures 		<- dataOptionsStructure$chosenClinicalFeatures
		chosenExtraClinicalFeatures <- dataOptionsStructure$chosenExtraClinicalFeatures
		chosenExpressionFeatures 	<- dataOptionsStructure$chosenExpressionFeatures
		extractionOutput 			<- DataExtractionCuratedOvarian(dataSet,chosenClinicalFeatures,chosenExtraClinicalFeatures,chosenExpressionFeatures,geneExpressionFlag,clinicalFeaturesFlag)
		allData 					<- extractionOutput$inputData
		allData 					<- allData[complete.cases(allData),,drop=FALSE]
		if(length(which(allData$days_to_event==0))!=0) allData <- allData[-which(allData$days_to_event==0),,drop=FALSE]
		modelVariables 				<- extractionOutput$allFeatureNames
		# Select random subset for bootstrapping
		allData 					<- allData[sample(nrow(allData), nrow(allData)),,drop=FALSE] 	# Randomise samples (to ensure folds are random)
		allX 						<- allData[,modelVariables,drop=FALSE]
		allY 						<- allData[,'days_to_event',drop=FALSE]
		allEvents 					<- allData[,'event',drop=FALSE]
		nSamples 					<- dim(allX)[1]													# Total number of samples
		dimension 					<- dim(allX)[2]
	} else if(dataSource=='TCGA2STAT'){
		dataSet 					<- dataOptionsStructure$dataSet
		geneExpressionFlag 			<- dataOptionsStructure$geneExpressionFlag
		clinicalFeaturesFlag 		<- dataOptionsStructure$clinicalFeaturesFlag
		chosenClinicalFeatures 		<- dataOptionsStructure$chosenClinicalFeatures
		chosenExpressionFeatures 	<- dataOptionsStructure$chosenExpressionFeatures
		nTraining 					<- dataOptionsStructure$nTraining
		nTest 						<- dataOptionsStructure$nTest
		extractionOutput 			<- DataExtractionTCGA2STAT(dataSet,chosenClinicalFeatures,chosenExpressionFeatures,clinicalFeaturesFlag,geneExpressionFlag)
		allData 					<- extractionOutput$inputData
		modelVariables 				<- extractionOutput$allFeatureNames
		# Deal with missing data
		onlyCompleteCases 	<- FALSE
		imputeMissingData 	<- TRUE
		imputeType 			<- 'ImputeBioc'
		if(onlyCompleteCases){
			allData 									<- allData[complete.cases(allData),,drop=FALSE]
		} else if(imputeMissingData){
			imputedStructure 							<- ImputeMissingData(allData[,modelVariables],imputeType)
			allData[,extractionOutput$allFeatureNames] 	<- imputedStructure$allData
		}
		if(sum(is.na(allData$days_to_event))!=0){
			allData 										<- allData[!is.na(allData$days_to_event),,drop=FALSE]
			if(nTraining+nTest>dim(allData)[1]) nTraining 	<-  dim(allData)[1] - nTest
		}
		if(sum(allData$days_to_event<=0)!=0){
			allData 										<- allData[!(allData$days_to_event<=0),,drop=FALSE]
			if(nTraining+nTest>dim(allData)[1]) nTraining 	<-  dim(allData)[1] - nTest
		}
		# Select random subset for bootstrapping
		allData 					<- allData[sample(nrow(allData), nrow(allData)),,drop=FALSE] 	# Randomise samples (to ensure folds are random)
		allX 						<- allData[,modelVariables,drop=FALSE]
		allY 						<- allData[,'days_to_event',drop=FALSE]
		allEvents 					<- allData[,'event',drop=FALSE]
		nSamples 					<- dim(allX)[1]													# Total number of samples
		dimension 					<- dim(allX)[2]
	} else if(dataSource=='TCGASynapse'){
		dataSet 					<- dataOptionsStructure$dataSet
		geneExpressionFlag 			<- dataOptionsStructure$geneExpressionFlag
		clinicalFeaturesFlag 		<- dataOptionsStructure$clinicalFeaturesFlag
		chosenClinicalFeatures 		<- dataOptionsStructure$chosenClinicalFeatures
		chosenExpressionFeatures 	<- dataOptionsStructure$chosenExpressionFeatures
		nTraining 					<- dataOptionsStructure$nTraining
		nTest 						<- dataOptionsStructure$nTest
		extractionOutput 			<- DataExtractionTCGASynapse(dataSet,chosenClinicalFeatures,chosenExpressionFeatures,clinicalFeaturesFlag,geneExpressionFlag)
		allData 					<- extractionOutput$inputData
		modelVariables 				<- extractionOutput$allFeatureNames
		# Deal with missing data
		onlyCompleteCases 	<- FALSE
		imputeMissingData 	<- TRUE
		imputeType 			<- 'ImputeBioc'
		if(onlyCompleteCases){
			allData 									<- allData[complete.cases(allData),,drop=FALSE]
		} else if(imputeMissingData){
			imputedStructure 							<- ImputeMissingData(allData[,modelVariables],imputeType)
			allData[,extractionOutput$allFeatureNames] 	<- imputedStructure$allData
		}
		if(sum(is.na(allData$days_to_event))!=0|sum(is.na(allData$event))!=0){
			allData 										<- allData[!is.na(allData$days_to_event)&!is.na(allData$event),,drop=FALSE]
			if(nTraining+nTest>dim(allData)[1]) nTraining 	<-  dim(allData)[1] - nTest
		}
		if(sum(allData$days_to_event<=0)!=0){
			allData 										<- allData[!(allData$days_to_event<=0),,drop=FALSE]
			if(nTraining+nTest>dim(allData)[1]) nTraining 	<-  dim(allData)[1] - nTest
		}
		# Select random subset for bootstrapping
		allData 					<- allData[sample(nrow(allData), nrow(allData)),,drop=FALSE] 	# Randomise samples (to ensure folds are random)
		allX 						<- allData[,modelVariables,drop=FALSE]
		allY 						<- allData[,'days_to_event',drop=FALSE]
		allEvents 					<- allData[,'event',drop=FALSE]
		nSamples 					<- dim(allX)[1]													# Total number of samples
		dimension 					<- dim(allX)[2]
	}

	trainingTestStructureNReps <- list()
	if(dataSource=='Generate'){
		for(i in 1:nReps){
			trainingIndices 			= sample(1:dataOptionsStructure$nSamples,nTraining)
			testIndices 				= sample((1:dataOptionsStructure$nSamples)[-trainingIndices],nTest)
			trainingData 				= as.matrix(allX[[i]][trainingIndices,,drop=FALSE])
			trainingTargets 			= as.matrix(allY[[i]][trainingIndices,,drop=FALSE])
			trainingEvents 				= as.matrix(allEvents[[i]][trainingIndices])
			testData 					= as.matrix(allX[[i]][testIndices,,drop=FALSE])
			testTargets 				= as.matrix(allY[[i]][testIndices,,drop=FALSE])
			testEvents 					= as.matrix(allEvents[[i]][testIndices])
			if(censoringType!='None'){
				trainingTargetsPreCensoring = allYPreCensoring[[i]][trainingIndices,,drop=FALSE]
				testTargetsPreCensoring 	= allYPreCensoring[[i]][testIndices,,drop=FALSE]
			}
			trainingTestStructureNReps[[i]] 	<- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,'testData'=testData,'testTargets'=testTargets,
														'events'=trainingEvents,'testEvents'=testEvents,'nTraining'=nTraining,'nTest'=nTest,'dimension'=dimension,
														'total'=dataOptionsStructure$nSamples,'methodFlag'=dataStructure[[1]]$methodFlag,
														'trainingTargetsPreCensoring'=trainingTargetsPreCensoring,'testTargetsPreCensoring'=testTargetsPreCensoring)
			if(censoringType!='None'){
				trainingTestStructureNReps[[i]]$trainingTargetsPreCensoring 	<- trainingTargetsPreCensoring
				trainingTestStructureNReps[[i]]$testTargetsPreCensoring 		<- testTargetsPreCensoring
			}
		}
	} else {
		for(i in 1:nReps){
			trainingIndices 			= sample(1:nSamples,nTraining)
			testIndices 				= sample((1:nSamples)[-trainingIndices],nTest)
			trainingData 				= as.matrix(allX[trainingIndices,,drop=FALSE])
			trainingTargets 			= as.matrix(allY[trainingIndices,,drop=FALSE])
			trainingEvents 				= as.matrix(allEvents[trainingIndices,,drop=FALSE])
			testData 					= as.matrix(allX[testIndices,,drop=FALSE])
			testTargets 				= as.matrix(allY[testIndices,,drop=FALSE])
			testEvents 					= as.matrix(allEvents[testIndices,,drop=FALSE])
			if(dataSource=='Generate'&censoringType!='None'){
				trainingTargetsPreCensoring = allYPreCensoring[trainingIndices,,drop=FALSE]
				testTargetsPreCensoring 	= allYPreCensoring[testIndices,,drop=FALSE]
			}
			trainingTestStructureNReps[[i]] 	<- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,'testData'=testData,'testTargets'=testTargets,
														'events'=trainingEvents,'testEvents'=testEvents,'nTraining'=nTraining,'nTest'=nTest,'dimension'=dimension,
														'total'=nSamples)
		}
	}

	return(trainingTestStructureNReps)
}

