GPRegression <- function(logHyp,parameterStructure,trainingTestStructure,logHypNoiseCorrection){
	#-------------------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------------------#

	#------------------------------------------- Extract data and parameters -------------------------------------------#
	trainingData 	<- trainingTestStructure$trainingData
	trainingTargets <- trainingTestStructure$trainingTargets
	trainingEvents 	<- trainingTestStructure$events
	nTraining 		<- trainingTestStructure$nTraining
	dimension 		<- trainingTestStructure$dimension

	noiseCorr 		<- parameterStructure$noiseCorr
	modelType 		<- parameterStructure$modelType
	optimType 		<- parameterStructure$optimType
	if(modelType=='survival') maxit <- parameterStructure$maxitSurvival else maxit <- parameterStructure$maxit
	meanFuncForm 	<- parameterStructure$meanFuncForm
	covFuncForm 	<- parameterStructure$covFuncForm
	maternParam 	<- parameterStructure$maternParam
	imposePriors 	<- parameterStructure$imposePriors

	logHypVec 		<- as.matrix(unlist(logHyp))
	if(meanFuncForm=='Zero') logHypVec <- logHypVec[-c((length(logHypVec)-dimension):length(logHypVec)),,drop=FALSE]
	if(noiseCorr=='noiseCorrLearned') logHypVec <- rbind(logHypVec,logHypNoiseCorrection)

	#---------------------------------------------- Learn hyperparameters ----------------------------------------------#
 	optimOutput 	<- optim(logHypVec,fn=ForOptimisationLog,gr=OptimisationGradientLog,trainingData=trainingData,trainingTargets=trainingTargets,
 							trainingEvents=trainingEvents,meanFuncForm=meanFuncForm,dimension=dimension,maternParam=maternParam,covFuncForm=covFuncForm,
 							nSamples=nTraining,noiseCorr=noiseCorr,logHypNoiseCorrection=logHypNoiseCorrection,imposePriors=imposePriors,
 							method=optimType,control=list('maxit'=maxit))
	logHypVec 		<- optimOutput$par
	objective 		<- optimOutput$value
	if(noiseCorr=='noiseCorrLearned'){
		logHypNoiseCorrection 	<- logHypVec[length(logHypVec),,drop=FALSE]
		logHypVec 				<- logHypVec[-length(logHypVec),,drop=FALSE]
	}
	varNoiseSqHyp 	<- exp(logHypVec[1,,drop=FALSE])
	varFuncSqHyp 	<- exp(logHypVec[2,,drop=FALSE])
	if(covFuncForm=='ARD'){
		lengthHyp 	<- exp(logHypVec[3:(2+dimension),,drop=FALSE])
	} else {
		lengthHyp 	<- exp(logHypVec[3,,drop=FALSE])
	}

	switch(meanFuncForm,
	    'Linear' = {meanHyp <- logHypVec[(length(logHypVec)-dimension):length(logHypVec),,drop=FALSE]}, 
	    'Zero'   = {meanHyp <- matrix(rep(0,dimension+1),nrow=dimension+1)}
	)

	#-------------- Compute mean vector, covariance matrix and log marginal likelihood using training data -------------#
	meanTraining 	<- MeanFunc(meanHyp,trainingData,meanFuncForm,nTraining)
	K 				<- CovFunc(t(trainingData),t(trainingData),maternParam,varFuncSqHyp,lengthHyp,covFuncForm)
	if(noiseCorr=='noiseCorrVec'|noiseCorr=='noiseCorrLearned'){
		varNoiseSqCorrection 					<- rep(0,nTraining)
		varNoiseSqCorrection[!trainingEvents] 	<- exp(logHypNoiseCorrection)
		L 										<- chol(K+diag(as.numeric(varNoiseSqHyp),dim(trainingData)[1],dim(trainingData)[1])+diag(as.numeric(varNoiseSqCorrection),dim(K)[1],dim(K)[2]))
	} else { 
		L     									<- chol(K+diag(as.numeric(varNoiseSqHyp),dim(K)[1],dim(K)[2]))
	}
	L 				<- t(L) 															# Make lower triangular, not upper triangular
	alpha 			<- solve(t(L),(solve(L,(trainingTargets-meanTraining))))
	logMarLik 		<- -(1/2)*t(trainingTargets-meanTraining)%*%alpha-sum(log(diag(L)))-(nTraining/2)*log(2*pi)

	logHypChosen 	<- list('noise'=log(varNoiseSqHyp),'func'=log(varFuncSqHyp),'length'=log(lengthHyp),'mean'=meanHyp)

	toReturn 		<- list('logHyp'=logHyp,'meanTraining'=meanTraining,'K'=K,'L'=L,'alpha'=alpha,'logHypChosen'=logHypChosen,'logMarLik'=logMarLik,'objective'=objective)
	if(noiseCorr=='noiseCorrVec'|noiseCorr=='noiseCorrLearned'){
		if(logHypNoiseCorrection==-Inf){
			logHypNoiseCorrection <- -10^10
		}
		toReturn$logHypNoiseCorrection <- logHypNoiseCorrection
	}
	return(toReturn)
}