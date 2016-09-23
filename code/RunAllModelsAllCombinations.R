RunAllModelsAllCombinations <- function(cancer,molPlatform,clinicalFlag,nReps){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	# Recreating analysis of TCGA data by Yuan et al., 2014
	# Cox PH, RSF and feature selection written using code from Synapse, ID:syn1720423, as reference. Where possible, original code has been used.
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# RunAllModelsAllCombinations applies each of 6 models to the provided data sets.
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	allParameterStructures              <- SetParametersRealYuan(cancer=cancer,molPlatform=molPlatform,clinicalFlag=clinicalFlag,nReps=nReps)
	parameterStructure                  <- allParameterStructures$parameterStructure
	plotSaveOptions                     <- allParameterStructures$plotSaveOptions
	dataOptionsStructure                <- allParameterStructures$dataOptionsStructure
	dataOptionsStructure                <- GetSynapseIDs(dataOptionsStructure)

	## Generate Data ##
	trainingTestStructureNReps    		<- GenerateDataYuanEtAl(dataOptionsStructure)

	parameterStructure$tolerance 		<- 0.00001*trainingTestStructureNReps[[1]]$nTraining*dataOptionsStructure$censoredProportion

	outputStructureCox  <- list()
	outputStructureRSF  <- list()
	outputStructureGPS1 <- list()
	outputStructureGP 	<- list()
	outputStructureGPS2 <- list()
	outputStructureGPS3 <- list()
	outputStructureGPS4 <- list()

	#####################################
	## Cox Proportional Hazards Model
	#####################################
	sink(file=paste0('Runs/',parameterStructure$unid,'/','output_',dataOptionsStructure$cancer,'_',dataOptionsStructure$molPlatform,'_',dataOptionsStructure$clinicalFlag,'_Cox','.txt'))
	for (seed in 1:dataOptionsStructure$nReps){
	    cat("seed =", seed,fill=TRUE)
	    cat("LASSO + cox: ",fill=TRUE)
	    outputStructureCox[[seed]]  <- ApplyCoxYuanEtAl(trainingTestStructureNReps[[seed]])
	    cat('-------------------------------------------------',fill=TRUE)
	}
	sink()

	#####################################
	## Random Survival Forest Model
	#####################################
	sink(file=paste0('Runs/',parameterStructure$unid,'/','output_',dataOptionsStructure$cancer,'_',dataOptionsStructure$molPlatform,'_',dataOptionsStructure$clinicalFlag,'_RSF','.txt'))
	for(seed in 1:dataOptionsStructure$nReps){
	    cat("seed =", seed,fill=TRUE)
	    cat("Random survival forest: ",fill=TRUE)
	    outputStructureRSF[[seed]]  <- ApplyRSFYuanEtAl(trainingTestStructureNReps[[seed]])
	    cat('-------------------------------------------------',fill=TRUE)
	}
	sink()

	#####################################
	## GP Without Censored Samples
	#####################################
	dataOptionsStructure$censoringType 		<- 'None'
	parameterStructure$noiseCorr 			<- FALSE
	parameterStructure$modelType 			<- 'non-survival'
	for(seed in 1:dataOptionsStructure$nReps){
	    trainingTestStructure 				<- trainingTestStructureNReps[[seed]]
		trainingTestStructure 				<- ApplyFeatureSelection(trainingTestStructure,dataOptionsStructure,parameterStructure)
		parameterStructure$logHypStart 		<- trainingTestStructure$logHypStart
		trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	    outputStructureGP[[seed]]  			<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	}

	#####################################
	## GP for Survival No Correction
	#####################################
	dataOptionsStructure$censoringType 		<- 'NormalLoopSample'
	parameterStructure$noiseCorr 			<- FALSE
	parameterStructure$modelType 			<- 'survival'
	for(seed in 1:dataOptionsStructure$nReps){
	    trainingTestStructure 				<- trainingTestStructureNReps[[seed]]
		trainingTestStructure 				<- ApplyFeatureSelection(trainingTestStructure,dataOptionsStructure,parameterStructure)
		parameterStructure$logHypStart 		<- trainingTestStructure$logHypStart
	    outputStructureGPS1[[seed]]  		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	}

	#####################################
	## GP for Survival V Correction
	#####################################
	dataOptionsStructure$censoringType 		<- 'NormalLoopSample'
	parameterStructure$noiseCorr 			<- 'noiseCorrVec'
	parameterStructure$modelType 			<- 'survival'
	for(seed in 1:dataOptionsStructure$nReps){
	    trainingTestStructure 				<- trainingTestStructureNReps[[seed]]
		trainingTestStructure 				<- ApplyFeatureSelection(trainingTestStructure,dataOptionsStructure,parameterStructure)
		parameterStructure$logHypStart 		<- trainingTestStructure$logHypStart
	    outputStructureGPS2[[seed]]  		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	}

	#####################################
	## GP for Survival L Correction
	#####################################
	dataOptionsStructure$censoringType 		<- 'NormalLoopSample'
	parameterStructure$noiseCorr 			<- 'noiseCorrLearned'
	parameterStructure$modelType 			<- 'survival'
	for(seed in 1:dataOptionsStructure$nReps){
	    trainingTestStructure 				<- trainingTestStructureNReps[[seed]]
		trainingTestStructure 				<- ApplyFeatureSelection(trainingTestStructure,dataOptionsStructure,parameterStructure)
		parameterStructure$logHypStart 		<- trainingTestStructure$logHypStart
	    outputStructureGPS3[[seed]]  		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	}

	toReturn <- list('outputStructureCox'=outputStructureCox,'outputStructureRSF'=outputStructureRSF,'outputStructureGP'=outputStructureGP,
					'outputStructureGPS1'=outputStructureGPS1,'outputStructureGPS2'=outputStructureGPS2,'outputStructureGPS3'=outputStructureGPS3,
					'dataOptionsStructure'=dataOptionsStructure,'parameterStructure'=parameterStructure)
	return(toReturn)
}
