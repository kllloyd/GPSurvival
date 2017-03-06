GPPredict <- function(regressionStructure,parameterStructure,trainingTestStructure,logHypNoiseCorrection){
	#-------------------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------------------#

	#------------------------------------------- Extract data and parameters -------------------------------------------#
	meanFuncForm 	<- parameterStructure$meanFuncForm
	covFuncForm 	<- parameterStructure$covFuncForm
	extraParam 		<- parameterStructure$extraParam
	noiseCorr 		<- parameterStructure$noiseCorr

	L 				<- regressionStructure$L
	alpha 			<- regressionStructure$alpha
	logHyp 			<- regressionStructure$logHypChosen

	trainingData 	<- trainingTestStructure$trainingData
	trainingTargets <- trainingTestStructure$trainingTargets
	trainingEvents 	<- trainingTestStructure$events
	nTraining 		<- trainingTestStructure$nTraining
	testData 		<- trainingTestStructure$testData
	nTest 			<- trainingTestStructure$nTest

	#--------------------------------------------- Convert log(hyp) to hyp ---------------------------------------------#
	varNoiseSqHyp 		<- exp(logHyp$noise)
	meanHyp 			<- logHyp$mean
	if(covFuncForm=='ARD'|covFuncForm=='SqExp'|covFuncForm=='Matern'|covFuncForm=='InformedARD'){
		varFuncSqHyp 	<- exp(logHyp$func)
		lengthHyp 		<- exp(logHyp$length)
	} else if(covFuncForm=='RF'){
		varFuncSqHyp 	<- NA
		lengthHyp 		<- NA
	}


	#---------------------------------------------- If needed, recompute L ---------------------------------------------#
	if(noiseCorr=='noiseCorrVec'|noiseCorr=='noiseCorrLearned'){
		varNoiseSqCorrection 					<- rep(0,nTraining)
		varNoiseSqCorrection[!trainingEvents] 	<- exp(logHypNoiseCorrection)
		trainingTargets 						<- trainingTestStructure$trainingTargets
		meanTraining 							<- MeanFunc(meanHyp,trainingData,meanFuncForm,nTraining)
		K 										<- CovFunc(trainingData,trainingData,trainingData,trainingTargets,extraParam,varFuncSqHyp,lengthHyp,covFuncForm)+diag(as.numeric(varNoiseSqHyp),dim(trainingData)[1],dim(trainingData)[1])
		K_noiseCorrected 						<- K+diag(as.numeric(varNoiseSqCorrection),dim(K)[1],dim(K)[2])
		L     									<- chol(K_noiseCorrected)
		L <- t(L)
		alpha 									<- solve(t(L),(solve(L,(trainingTargets-meanTraining))))
	}

	#---------------------------- Compute mean vector and covariance matrix using test data ----------------------------#
	meanTest 		<- MeanFunc(meanHyp,testData,meanFuncForm,nTest)
	kStar 			<- CovFunc(trainingData,testData,trainingData,trainingTargets,extraParam,varFuncSqHyp,lengthHyp,covFuncForm)
	v 				<- solve(L,kStar)

	#------- Predict function mean, function variance and data variance (data mean = function mean) at test point ------#
	funcMeanPred 	<- meanTest+t(kStar)%*%alpha
	varFunction  	<- diag(CovFunc(testData,testData,trainingData,trainingTargets,extraParam,varFuncSqHyp,lengthHyp,covFuncForm)-t(v)%*%v)
	varData      	<- varFunction+as.numeric(varNoiseSqHyp)

	toReturn 		<- list('funcMeanPred'=funcMeanPred,'varFunction'=varFunction,'varData'=varData)
	return(toReturn)
}