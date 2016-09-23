PreLearnHyperparam <- function(parameterStructure,plotSaveOptions,preLearnStructure){
    #-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Function pre-learns hyperparameter values using subset of non-censored training samples.
	# Fast burn-in of hyperparameter values to save time when updating censored training targets
	#-------------------------------------------------------------------------------------------------------#

	source('ForOptimisationLog.R')
	source('OptimisationGradientLog.R')

	meanFuncForm 	<- parameterStructure$meanFuncForm
	maternParam 	<- parameterStructure$maternParam
	covFuncForm 	<- parameterStructure$covFuncForm
	logHypStart 	<- parameterStructure$logHypStart
	optimType 		<- parameterStructure$optimType
	maxitPreLearn 	<- parameterStructure$maxitPreLearn
	noiseCorr 		<- parameterStructure$noiseCorr
	model 			<- parameterStructure$model

	saveFiles 		<- plotSaveOptions$saveFiles
	printResults 	<- plotSaveOptions$printResults
	fileName 		<- plotSaveOptions$fileName
	savePlots 		<- plotSaveOptions$savePlots
	printPlots 		<- plotSaveOptions$printPlots

	trainingData 	<- preLearnStructure$trainingData
	trainingTargets <- preLearnStructure$trainingTargets
	trainingEvents 	<- preLearnStructure$trainingEvents
	testData 		<- preLearnStructure$testData
	testTargets 	<- preLearnStructure$testTargets
	dimension 		<- preLearnStructure$dimension
	nTraining 		<- preLearnStructure$nTraining
	nTest 			<- preLearnStructure$nTest
	applyCorr		<- preLearnStructure$applyCorr

	#---------------------------- Extract hyperparameters ----------------------------#
	logHypVec 		<- as.matrix(unlist(logHypStart))
	if(meanFuncForm=='Zero') logHypVec <- logHypVec[-c((length(logHypVec)-dimension):length(logHypVec)),,drop=FALSE]

	#----------------------------- Learn hyperparameters -----------------------------#
	optimOutput 	<- optim(logHypVec,fn=ForOptimisationLog,gr=OptimisationGradientLog,trainingData=trainingData,trainingTargets=trainingTargets,
							trainingEvents=trainingEvents,meanFuncForm=meanFuncForm,dimension=dimension,maternParam=maternParam,covFuncForm=covFuncForm,
	                        nSamples=nTraining,noiseCorr=FALSE,logHypNoiseCorrection=NA,method=optimType,control=list('maxit'=maxitPreLearn)) #,'trace'=9
	logHypVec 		<- optimOutput$par
	logHypNoise 	<- logHypVec[1,,drop=FALSE]
	logHypFunc 		<- logHypVec[2,,drop=FALSE]
	if(covFuncForm=='ARD'){
		logHypLength 		<- logHypVec[3:(2+dimension),,drop=FALSE]
	} else {logHypLength 	<- logHypVec[3,,drop=FALSE]}
	switch(meanFuncForm,
	    'Linear' = {logHypMean 	<- logHypVec[(length(logHypVec)-dimension):length(logHypVec),,drop=FALSE]}, 
	    'Zero'  = {logHypMean 	<- matrix(rep(0,dimension+1),nrow=dimension+1)}
	)
	logHypLearned <- list('noise'=logHypNoise,'func'=logHypFunc,'length'=logHypLength,'mean'=logHypMean)
	# cat(unlist(logHypLearned),fill=TRUE)
	if(model=='GPSurvLaplace'){
		learningDataStructure 	<- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,
										'testData'=testData,'nTraining'=nTraining,'nTest'=nTest,
										'dimension'=dimension)
		regressionStructure 	<- GPRegression(logHypLearned,parameterStructure,learningDataStructure,NA)
		finalDataStructure 		<- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,'events'=trainingEvents,
										'testData'=testData,'nTraining'=nTraining,'nTest'=nTest,
										'dimension'=dimension)
		predictionStructure 	<- GPPredict(regressionStructure,parameterStructure,finalDataStructure,NA)
		fhatBurnIn 				<- predictionStructure$funcMeanPred
		fhatBurnInVarFunction  	<- predictionStructure$varFunction
		fhatBurnInVarData      	<- predictionStructure$varData
		cat(fhatBurnIn,fill=TRUE)
	}

	#---------- Print and save pre- and post-learning hyperparameter values ----------#
	if(printResults){
		cat('-------------------------------------------------------',fill=TRUE)
		cat('Hyperparameter initialisation values:',fill=TRUE)
		cat('\tlog(noise hyp) =',logHypStart$noise,fill=TRUE)
		cat('\tlog(function hyp) =',logHypStart$func,fill=TRUE)
		cat('\tlog(length hyp) =',logHypStart$length,fill=TRUE)
		cat('\tmean hyp =',if(meanFuncForm=='Zero') rep(0,dimension+1) else logHypStart$mean,fill=TRUE)
		cat('Hyperparameter pre-learned values:',fill=TRUE)
		cat('\tlog(noise hyp) =',logHypLearned$noise,fill=TRUE)
		cat('\tlog(function hyp) =',logHypLearned$func,fill=TRUE)
		cat('\tlog(length hyp) =',logHypLearned$length,fill=TRUE)
		cat('\tmean hyp =',if(meanFuncForm=='Zero') rep(0,dimension+1) else logHypLearned$mean,fill=TRUE)
	}

	if(saveFiles){
		sink(file=paste0(fileName,'BurnInHyp.txt'))
			cat('-------------------------------------------------------',fill=TRUE)
			cat('Hyperparameter initialisation values:',fill=TRUE)
			cat('\tlog(noise hyp) =',logHypStart$noise,fill=TRUE)
			cat('\tlog(function hyp) =',logHypStart$func,fill=TRUE)
			cat('\tlog(length hyp) =',logHypStart$length,fill=TRUE)
			cat('\tmean hyp =',if(meanFuncForm=='Zero') rep(0,dimension+1) else logHypStart$mean,fill=TRUE)
			cat('Hyperparameter pre-learned values:',fill=TRUE)
			cat('\tlog(noise hyp) =',logHypLearned$noise,fill=TRUE)
			cat('\tlog(function hyp) =',logHypLearned$func,fill=TRUE)
			cat('\tlog(length hyp) =',logHypLearned$length,fill=TRUE)
			cat('\tmean hyp =',if(meanFuncForm=='Zero') rep(0,dimension+1) else logHypLearned$mean,fill=TRUE)
		sink()
	}

	if((printPlots|savePlots)&model=='GPSurvLaplace'){
		if(dimension==1){
			xLine 					<- matrix(seq(from=min(c(trainingData,testData)),to=max(c(trainingData,testData)),length.out=100),ncol=1)
			finalDataStructure2 	<- list('trainingData'=trainingData,'trainingTargets'=trainingTargets,
											'nTraining'=nTraining,'testData'=xLine,'nTest'=100,'dimension'=dimension,
											'events'=trainingEvents)
			predictionStructure2 	<- GPPredict(regressionStructure,parameterStructure,finalDataStructure2,NA)
			meanLine 				<- predictionStructure2$funcMeanPred
			varLine      			<- predictionStructure2$varData
			f 						<- rbind(as.matrix(meanLine+2*sqrt(varLine)),as.matrix(rev(meanLine-2*sqrt(varLine))))
			plot(c(0,0),ylim=c(min(c(testTargets,trainingTargets,fhatBurnIn)),max(c(testTargets,trainingTargets,fhatBurnIn))),xlim=c(min(testData,trainingData),max(testData,trainingData)),col='white',ylab='Targets',xlab='Data',main=paste0('hypChosen = [',round(exp(logHypLearned$noise),4),',',round(exp(logHypLearned$func),4),',',round(exp(logHypLearned$length),4),']'))
			polygon(rbind(xLine,as.matrix(rev(xLine))),f,col=rgb(173/250,216/250,230/250,0.5))
			lines(xLine,meanLine)
			points(testData,testTargets,pch=20,col='black')
			points(testData,fhatBurnIn,pch=20,col='blue')
			points(trainingData,trainingTargets,pch=20,col='green')
			segments(x0=testData,y0=testTargets,x1=testData,y1=fhatBurnIn,col='black')
			plotPreLearn <- recordPlot()
		}
	} else {
		plotPreLearn <- NA
	}

	toReturn <- list('logHypLearned'=logHypLearned)
	if(model=='GPSurvLaplace'){
		toReturn$fhatBurnIn 			<- fhatBurnIn
		toReturn$fhatBurnInVarFunction 	<- fhatBurnInVarFunction
		toReturn$fhatBurnInVarData 		<- fhatBurnInVarData
		if((printPlots|savePlots)&dimension==1) toReturn$plotPreLearn <- plotPreLearn
	}

	return(toReturn)
}