ApplyGPDiffInf <- function(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Hyperparameter values may be pre-learned using a subset of un-censored samples (burn-in).
	# A model is trained on the whole training set and used to predict censored training set survival times. 
	# Training set targets are updated. This is repeated until the process converges.
	# Predictions are then made on the test set.
	# Gaussian Process: Options for covariance function and mean function, and hyperparameter starting values
	#                   For covariance function only 'SqExp' is currently working ('Matern' will also be implemented and 'ARD' is in progress)
	#                   Mean function options are 'Zero' or 'Linear'. Linear uses more hyperparameters.
	# 					Noise variance correction may be none (GPSurvNoCorr = GPS1), 'noiseCorrLearned' (GPSurvCorrL = GPS2) or 'noiseCorrVec' (GPSurvCorrV = GPS3)
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------- Load required libraries and functions --------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	source('ForOptimisationLog.R')
	source('OptimisationGradientLog.R')
	source('GPRegression.R')
	source('GPPredict.R')
	source('CovFunc.R')
	source('MeanFunc.R')
	source('AdjustTrainingSurvivalMeanVariance.R')

	modelType 		<- parameterStructure$modelType
	tolerance 		<- parameterStructure$tolerance
	maxCount 		<- parameterStructure$maxCount
	hypChangeTol 	<- parameterStructure$hypChangeTol
	unid 			<- parameterStructure$unid
	noiseCorr 		<- parameterStructure$noiseCorr
	inferenceType 	<- parameterStructure$inferenceType
	burnIn 			<- parameterStructure$burnIn
	logHypStart 	<- parameterStructure$logHypStart
	outerFolder 	<- plotSaveOptions$outerFolder
	censoringType 	<- dataOptionsStructure$censoringType

	if(modelType=='survival'){
		if(inferenceType=='median'){
			model <- 'GPSurvInfMed'
		} else if(inferenceType=='uniform'){
			model <- 'GPSurvInfUnif'
		}
	} else {
			model <- 'GPNonSurv'
	}
	parameterStructure$model <- model

	#-------------------------------------------------------------------------------------------------------#
	#--------------------------- Flags for printing and saving options and plots ---------------------------#
	#-------------------------------------------------------------------------------------------------------#
	plotRuns 		<- plotSaveOptions$plotRuns
	printResults 	<- plotSaveOptions$printResults
	savePlots 		<- plotSaveOptions$savePlots
	saveFiles 		<- plotSaveOptions$saveFiles
	printPlots 		<- plotSaveOptions$printPlots

	fileName 		<- paste0(getwd(),'/',outerFolder,'/',unid,'/',model,'/',unid,model)

	if(savePlots|saveFiles) {
		dir.create(file.path(paste0(getwd(),"/",outerFolder),unid),showWarnings=FALSE)
		dir.create(file.path(paste0(getwd(),"/",outerFolder,'/',unid),model),showWarnings=FALSE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#---------------------------------- Printing model options to screen -----------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	## Put code here to print model options to screen. Normally turned off using printOptions.
	# if(printOptions){
		
	# }


	#-------------------------------------------------------------------------------------------------------#
	#------------------------------ Separate censored and uncensored samples -------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	if(modelType=='survival'){
		trainingTestStructure$trainingTargetsCensored 		<- trainingTestStructure$trainingTargets[trainingTestStructure$events==0,,drop=FALSE]
		trainingTestStructure$trainingTargetsDied 			<- trainingTestStructure$trainingTargets[trainingTestStructure$events==1,,drop=FALSE]
		trainingTestStructure$trainingDataCensored 			<- trainingTestStructure$trainingData[trainingTestStructure$events==0,,drop=FALSE]
		trainingTestStructure$trainingDataDied 				<- trainingTestStructure$trainingData[trainingTestStructure$events==1,,drop=FALSE]
		trainingTestStructure$nTrainingCensored 			<- dim(trainingTestStructure$trainingDataCensored)[1]
		trainingTestStructure$nTrainingDied 				<- dim(trainingTestStructure$trainingDataDied)[1]
		trainingTestStructure$trainingTargetsCensoredStart 	<- trainingTestStructure$trainingTargetsCensored
	}

	timeStart <- Sys.time()


	#-------------------------------------------------------------------------------------------------------#
	#-------------- Use inference method to predict survival of censored training set members --------------#
	#-------------------------------------------------------------------------------------------------------#
	if(inferenceType=='uniform'){
		trainingTestStructure$trainingTargetsLearned 	<- trainingTestStructure$trainingTargets
		for(i in 1:dim(trainingTestStructure$trainingTargets)[1]){
			if(trainingTestStructure$events[i]==0) trainingTestStructure$trainingTargetsLearned[i,1] <- runif(1,min=trainingTestStructure$trainingTargets[i,1],max=max(trainingTestStructure$trainingTargets))
		}
	} else if(inferenceType=='median'){
		trainingTestStructure$trainingTargetsLearned 									<- trainingTestStructure$trainingTargets
		trainingTestStructure$trainingTargetsLearned[trainingTestStructure$events==0] 	<- median(trainingTestStructure$trainingTargetsDied)
	}


	#-------------------------------------------------------------------------------------------------------#
	#--------------------------- Train GP to predict survival of test set members --------------------------#
	#-------------------------------------------------------------------------------------------------------#
	learningDataStructure 			<- list('trainingData'=trainingTestStructure$trainingData,
											'trainingTargets'=trainingTestStructure$trainingTargetsLearned,
											'testData'=trainingTestStructure$testData,'nTraining'=trainingTestStructure$nTraining,
											'nTest'=trainingTestStructure$nTest,'dimension'=trainingTestStructure$dimension)
	parameterStructure$noiseCorr 	<- FALSE
	logHyp 							<- logHypStart
	regressionStructure 			<- GPRegression(logHyp,parameterStructure,learningDataStructure,NA)
	logHypNoiseCorrection 			<- NA
	predictionStructure 			<- GPPredict(regressionStructure,parameterStructure,trainingTestStructure,logHypNoiseCorrection)

	logHypChosen 					<- regressionStructure$logHypChosen
	funcMeanPred 					<- predictionStructure$funcMeanPred
	varFunction  					<- predictionStructure$varFunction
	varData      					<- predictionStructure$varData

	if((printPlots|savePlots)&trainingTestStructure$dimension==1){
		xLine 					<- matrix(seq(from=min(c(trainingTestStructure$trainingData,trainingTestStructure$testData)),to=max(c(trainingTestStructure$trainingData,trainingTestStructure$testData)),length.out=100),ncol=1)
		finalDataStructure2 	<- list('trainingData'=trainingTestStructure$trainingData,
										'trainingTargets'=trainingTestStructure$trainingTargetsLearned,
										'testData'=xLine,
										'nTraining'=trainingTestStructure$nTraining,
										'nTest'=100,'dimension'=trainingTestStructure$dimension,
										'events'=trainingTestStructure$events)
		predictionStructure2 	<- GPPredict(regressionStructure,parameterStructure,finalDataStructure2,logHypNoiseCorrection)
		meanLine 				<- predictionStructure2$funcMeanPred
		varLine      			<- predictionStructure2$varData
	}

	if(modelType=='survival'){
		## Predictions on training set ##
		trainPredDataStructure 	<- list('trainingData'=trainingTestStructure$trainingData,
									'trainingTargets'=trainingTestStructure$trainingTargetsLearned,
									'testData'=trainingTestStructure$trainingData,'nTraining'=trainingTestStructure$nTraining,
									'nTest'=trainingTestStructure$nTraining,'dimension'=trainingTestStructure$dimension,
									'events'=trainingTestStructure$events)
		predictionStructure3 	<- GPPredict(regressionStructure,parameterStructure,trainPredDataStructure,logHypNoiseCorrection)
		trainingSetPredictions 	<- predictionStructure3$funcMeanPred
	} else {
		trainPredDataStructure 	<- list('trainingData'=trainingTestStructure$trainingData,
									'trainingTargets'=trainingTestStructure$trainingTargets,
									'testData'=trainingTestStructure$trainingData,'nTraining'=trainingTestStructure$nTraining,
									'nTest'=trainingTestStructure$nTraining,'dimension'=trainingTestStructure$dimension,
									'events'=trainingTestStructure$events)
		predictionStructure3 	<- GPPredict(regressionStructure,parameterStructure,trainPredDataStructure,logHypNoiseCorrection)
		trainingSetPredictions 	<- predictionStructure3$funcMeanPred
	}

	#-------------------------------------------------------------------------------------------------------#
	#----------------------------------------- Performance Metrics -----------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	if(modelType=='survival'){
		if(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse'|dataOptionsStructure$dataSource=='TCGAYuanEtAl'){
			metricStructure <- CalculateMetrics(funcMeanPred,trainingTestStructure$testTargets,trainingTestStructure$testEvents,trainingSetPredictions,trainingTestStructure$trainingTargets,trainingTestStructure$events)
			c.index 		<- metricStructure$c.index
			rmse 			<- metricStructure$rmse
		} else {
			metricStructure <- CalculateMetrics(funcMeanPred,trainingTestStructure$testTargetsPreCensoring,rep(1,dim(funcMeanPred)[1]),trainingSetPredictions,trainingTestStructure$trainingTargets,trainingTestStructure$events)
			c.index 		<- metricStructure$c.index
			rmse 			<- metricStructure$rmse
		}
	} else {
		metricStructure <- CalculateMetrics(funcMeanPred,trainingTestStructure$testTargets,rep(1,dim(funcMeanPred)[1]),trainingSetPredictions,trainingTestStructure$trainingTargets,rep(1,dim(trainingSetPredictions)[1]))
		c.index 		<- metricStructure$c.index
		rmse 			<- metricStructure$rmse
	}

	timeEnd   			<- Sys.time()
	timeTaken 			<- difftime(timeEnd,timeStart, units='min')

	#-------------------------------------------------------------------------------------------------------#
	#--------------------------------- Print and save variables and plots ----------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	if(printResults){
		cat('-------------------------------------------------------',fill=TRUE)
		cat('Model =',model,fill=TRUE)
		cat('Time taken = ',timeTaken,fill=TRUE)
		cat('C Index =',round(c.index,4),fill=TRUE)
		cat('RMSE =',round(rmse,4),fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}

	if(saveFiles){
		sink(paste0(fileName,'Metrics.txt'),append=TRUE)
			cat('Model =',model,fill=TRUE)
			cat('Time taken = ',timeTaken,fill=TRUE)
			cat('RMSE =',rmse,fill=TRUE)
			cat('Concordance Index =',c.index,fill=TRUE)
			cat('-------------------------------------------------------',fill=TRUE)
		sink()

		sink(paste0(fileName,'LearnHyperparam.txt'),append=TRUE)
			cat('Hyperparameter initialisation values:',fill=TRUE)
			cat('\tlog(noise hyp) =',logHypStart$noise,fill=TRUE)
			cat('\tlog(function hyp) =',logHypStart$func,fill=TRUE)
			cat('\tlog(length hyp) =',logHypStart$length,fill=TRUE)
			cat('\tmean hyp =',if(parameterStructure$meanFuncForm=='Zero') rep(0,trainingTestStructure$dimension+1) else logHypStart$mean,fill=TRUE)
			cat('Hyperparameter final values:',fill=TRUE)
			cat('\tlog(noise hyp) =',regressionStructure$logHypChosen$noise,fill=TRUE)
			cat('\tlog(function hyp) =',regressionStructure$logHypChosen$func,fill=TRUE)
			cat('\tlog(length hyp) =',regressionStructure$logHypChosen$length,fill=TRUE)
			cat('\tmean hyp =',if(parameterStructure$meanFuncForm=='Zero') rep(0,trainingTestStructure$dimension+1) else regressionStructure$logHypChosen$mean,fill=TRUE)
			cat('-------------------------------------------------------',fill=TRUE)
		sink()

		write.table(funcMeanPred,paste0(fileName,'TestTargetPredictions.csv'),sep=',',quote=FALSE,row.names=FALSE,append=TRUE)
		if(censoringType!='None') write.table(trainingTestStructure$testTargetsPreCensoring,paste0(fileName,'testTargetsPreCensoring.csv'),sep=',',quote=FALSE,row.names=FALSE,append=TRUE)
	}

	if(printPlots|savePlots){
		if(trainingTestStructure$dimension==1){
			if(modelType=='survival'&censoringType!='None'){
				opar <- par('mar')
				layout(rbind(1,2), heights=c(7,1))
			}
			trainingData 		<- trainingTestStructure$trainingData
			trainingTargets 	<- trainingTestStructure$trainingTargets
			testData 			<- trainingTestStructure$testData
			f 					<- rbind(as.matrix(meanLine+2*sqrt(varLine)),as.matrix(rev(meanLine-2*sqrt(varLine))))
			if(modelType=='survival'&censoringType!='None'){
				plot(c(0,0),ylim=c(min(c(trainingTestStructure$trainingTargets,trainingTestStructure$trainingTargetsPreCensoring,trainingTargetsCensoredEachRun[,dim(trainingTargetsCensoredEachRun)[2]])),max(c(trainingTestStructure$trainingTargets,trainingTestStructure$trainingTargetsPreCensoring,trainingTargetsCensoredEachRun[,dim(trainingTargetsCensoredEachRun)[2]]))),xlim=c(min(trainingTestStructure$trainingData,trainingTestStructure$testData),max(trainingTestStructure$trainingData,trainingTestStructure$testData)),col='white',ylab='Targets',xlab='Data',main=paste0('hypChosen = [',round(exp(logHypChosen$noise),4),',',round(exp(logHypChosen$func),4),',',round(exp(logHypChosen$length),4),']'))
			} else {
				plot(c(0,0),ylim=c(min(trainingTestStructure$trainingTargets),max(trainingTestStructure$trainingTargets)),xlim=c(min(trainingTestStructure$trainingData,trainingTestStructure$testData),max(trainingTestStructure$trainingData,trainingTestStructure$testData)),col='white',ylab='Targets',xlab='Data',main=paste0('hypChosen = [',round(exp(logHypChosen$noise),4),',',round(exp(logHypChosen$func),4),',',round(exp(logHypChosen$length),4),']'))
			}
			polygon(rbind(xLine,as.matrix(rev(xLine))),f,col=rgb(173/250,216/250,230/250,0.5))
			lines(xLine,meanLine)
			points(trainingTestStructure$trainingData[trainingTestStructure$events==0],trainingTestStructure$trainingTargets[trainingTestStructure$events==0],pch=8,col='black')
			if(censoringType!='None') points(trainingTestStructure$trainingData,trainingTestStructure$trainingTargetsPreCensoring,pch=1,col='blue')
			if(censoringType=='None') points(trainingTestStructure$trainingData,trainingTestStructure$trainingTargets,pch=1,col='blue')
			if(modelType=='survival') points(trainingTestStructure$trainingData[trainingTestStructure$events==0],trainingTestStructure$trainingTargetsLearned[trainingTestStructure$events==0],pch=6,col='green')
			points(testData,funcMeanPred,col='purple',pch=18)
			if(modelType=='survival') segments(x0=trainingTestStructure$trainingData[trainingTestStructure$events==0],y0=trainingTestStructure$trainingTargets[trainingTestStructure$events==0],x1=trainingTestStructure$trainingData[trainingTestStructure$events==0],y1=trainingTestStructure$trainingTargetsLearned[trainingTestStructure$events==0],col='black')
			if(modelType=='survival'&censoringType!='None'){
				par(mar=c(0,0,0,0))
				plot.new()
				legend('center',c('Training, Pre-censoring','Training, Post-censoring','Training, Learned','Test'),col=c('blue','black','green','purple'),pch=c(1,8,6,18),lty=c(0,0,0,0),ncol=2,bty ="n")
			}
			plot1 <- recordPlot()
			if(modelType=='survival'&censoringType!='None'){
				par(mar=opar)
				layout(1)
			}
		} else if(trainingTestStructure$dimension==2){
			plot1 <- NULL
			# surface3d(interp(trainingTestStructure$testData[,1],trainingTestStructure$testData[,2],funcMeanPred),col=rainbow(100)[cut(funcMeanPred,100)],back='lines',front='lines')
			# plot3d(trainingTestStructure$trainingData[,1],trainingTestStructure$trainingData[,2],trainingTestStructure$trainingTargets,type='p',add=TRUE)
			# plot3d(trainingTestStructure$trainingData[,1],trainingTestStructure$trainingData[,2],trainingTestStructure$trainingTargetsPreCensoring,type='p',add=TRUE,col='blue')
		} else {
			plot1 <- NULL
		}

		if(censoringType!='None'&(dataOptionsStructure$dataSource!='CuratedOvarian'&dataOptionsStructure$dataSource!='TCGA2STAT'&dataOptionsStructure$dataSource!='TCGASynapse'&dataOptionsStructure$dataSource!='TCGAYuanEtAl')){
			plot(c(1,1),col='white',main=paste0('c index = ',round(c.index,4),', rmse = ',round(rmse,4)),xlab='Measured',ylab='Predicted',ylim=c(min(c(trainingTestStructure$testTargetsPreCensoring,funcMeanPred)),max(c(trainingTestStructure$testTargetsPreCensoring,funcMeanPred))),xlim=c(min(c(trainingTestStructure$testTargetsPreCensoring,funcMeanPred)),max(c(trainingTestStructure$testTargetsPreCensoring,funcMeanPred))))
			points(trainingTestStructure$testTargetsPreCensoring,funcMeanPred)
			abline(0,1)
			plot2 <- recordPlot()
		} else if(censoringType!='None'&(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse'|dataOptionsStructure$dataSource=='TCGAYuanEtAl')){
			plot(c(1,1),col='white',main=paste0('c index = ',round(c.index,4),', rmse = ',round(rmse,4)),xlab='Measured',ylab='Predicted',ylim=c(min(c(trainingTestStructure$testTargets,funcMeanPred)),max(c(trainingTestStructure$testTargets,funcMeanPred))),xlim=c(min(c(trainingTestStructure$testTargets,funcMeanPred)),max(c(trainingTestStructure$testTargets,funcMeanPred))))
			points(trainingTestStructure$testTargets[trainingTestStructure$testEvents==0],funcMeanPred[trainingTestStructure$testEvents==0],col='green')
			points(trainingTestStructure$testTargets[trainingTestStructure$testEvents==1],funcMeanPred[trainingTestStructure$testEvents==1],col='blue')
			abline(0,1)
			plot2 <- recordPlot()
		} else {
			plot(c(1,1),col='white',main=paste0('c index = ',round(c.index,4),', rmse = ',round(rmse,4)),xlab='Measured',ylab='Predicted',ylim=c(min(c(trainingTestStructure$testTargets,funcMeanPred)),max(c(trainingTestStructure$testTargets,funcMeanPred))),xlim=c(min(c(trainingTestStructure$testTargets,funcMeanPred)),max(c(trainingTestStructure$testTargets,funcMeanPred))))
			points(trainingTestStructure$testTargets,funcMeanPred)
			abline(0,1)
			plot2 <- recordPlot()
		}

		if(censoringType!='None'&modelType=='survival'){
			if(dataOptionsStructure$dataSource=='Generate'){
				PlotKaplanMeier(funcMeanPred,trainingTestStructure$testTargetsPreCensoring,rep(1,dim(trainingTestStructure$testTargets)[1]),model)
				plot3 <- recordPlot()
			} else {
				PlotKaplanMeier(funcMeanPred,trainingTestStructure$testTargets,trainingTestStructure$testEvents,model)
				plot3 <- recordPlot()
			}
		} else {
			plot3 <- NULL
		}
	}

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------------------- Return output --------------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	toReturn 														<- predictionStructure
	toReturn$parameterStructure 									<- parameterStructure
	toReturn$dataOptionsStructure 									<- dataOptionsStructure
	toReturn$trainingTestStructure		 							<- trainingTestStructure
	toReturn$trainingSetPredictions 								<- trainingSetPredictions
	toReturn$logHypChosen 											<- regressionStructure$logHypChosen
	toReturn$c.index 												<- c.index
	toReturn$rmse 													<- rmse
	toReturn$timeTaken												<- timeTaken
	if(printPlots|savePlots) toReturn$plot1 						<- plot1
	if(printPlots|savePlots) toReturn$plot2 						<- plot2
	if(printPlots|savePlots) toReturn$plot3 						<- plot3
	if((printPlots|savePlots)&trainingTestStructure$dimension==1){
		toReturn$modelPlotOutputStructure			<- modelPlotOutputStructure
		toReturn$meanLine 							<- meanLine
		toReturn$varLine 							<- varLine
		toReturn$xLine 								<- xLine
	}
	return(toReturn)

}
