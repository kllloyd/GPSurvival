RemoveCensored <- function(trainingTestStructure,dataOptionsStructure){
	#---------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#---------------------------------------------------------------------------------#
	# Removes censored samples from the data set
	#---------------------------------------------------------------------------------#

	dataSource 		<- dataOptionsStructure$dataSource

	trainingEvents 	<- trainingTestStructure$events
	testEvents 		<- trainingTestStructure$testEvents
	
	trainingData 	<- trainingTestStructure$trainingData
	trainingTargets <- trainingTestStructure$trainingTargets
	testData 		<- trainingTestStructure$testData
	testTargets 	<- trainingTestStructure$testTargets
	if(dataSource!='CuratedOvarian'&dataSource!='TCGA2STAT'&dataSource!='TCGASynapse'&dataSource!='TCGAYuanEtAl'){
		testTargetsPreCensoring <- trainingTestStructure$testTargetsPreCensoring
	}
	
	trainingData 	<- trainingData[trainingEvents==1,,drop=FALSE]
	trainingTargets <- trainingTargets[trainingEvents==1,,drop=FALSE]
	if(dataSource=='CuratedOvarian'|dataSource=='TCGA2STAT'|dataSource=='TCGASynapse'|dataSource=='TCGAYuanEtAl'){
		testData 		<- testData[testEvents==1,,drop=FALSE]
		testTargets 	<- testTargets[testEvents==1,,drop=FALSE]
	} else {
		testTargets 	<- testTargetsPreCensoring
	}
	
	nTraining 		<- dim(trainingData)[1]
	nTest 			<- dim(testData)[1]
	dimension 		<- dim(trainingData)[2]

	trainingEvents 	<- rep(1,nTraining)
	testEvents 		<- rep(1,nTest)

	trainingTestStructure <- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,'testData'=testData,'testTargets'=testTargets,'nTraining'=nTraining,
								'nTest'=nTest,'dimension'=dimension,'events'=trainingEvents,'testEvents'=testEvents)

	return(trainingTestStructure)
}