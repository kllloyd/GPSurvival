ApplyAFT <- function(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions){
	#-------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#-------------------------------------------------------------------------------------------------------#
	# Predicting risk of data set members given censored and uncensored training data.
	# Applied Accelerated Failure Time model using survreg
	#-------------------------------------------------------------------------------------------------------#

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------- Load required libraries and functions --------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	library('survival')
	library('MASS')
	library('glmnet')
	library('randomForestSRC')
	library('gbm')
	library('zoo')
	library('survcomp')
	library('ipred')

	#-------------------------------------------------------------------------------------------------------#
	#--------------------------- Flags for printing and saving options and plots ---------------------------#
	#-------------------------------------------------------------------------------------------------------#
	printOptions    = plotSaveOptions$printOptions
	printResults 	= plotSaveOptions$printResults
	printPlots 		= plotSaveOptions$printPlots
	saveFiles 		= plotSaveOptions$saveFiles
	savePlots 		= plotSaveOptions$savePlots

	model 			= 'AFT'
	unid 			= parameterStructure$unid
	folderName 		= dataOptionsStructure$folderName
	outerFolder 	= plotSaveOptions$outerFolder
	fileName 		= paste0(getwd(),'/',folderName,'/',model,'/',unid,model)
	if(savePlots|saveFiles) {
		dir.create(file.path(paste0(getwd(),"/",outerFolder),unid),showWarnings=FALSE)
		dir.create(file.path(paste0(getwd(),"/",outerFolder,'/',unid),model),showWarnings=FALSE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------------------- Model options --------------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	if(printOptions){
		cat('-------------------------------------------------------',fill=TRUE)
		cat('Predicting y from x data by generating synthetic data.',fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
		cat('Data generation options:',fill=TRUE)
		cat('\tGeneration was using Gaussian processes',fill=TRUE)
		cat(paste0('\t',dataOptionsStructure$nSamples),'samples were generated with',dataOptionsStructure$dimension,'dimensions',fill=TRUE)
		cat('\tGenerating log hyperparameters were',unlist(dataOptionsStructure$logHypGenerate),fill=TRUE)
		cat('\tSamples were produced on a grid between',dataOptionsStructure$gridMinimum,'and',dataOptionsStructure$gridMaximum,fill=TRUE)
		cat('\tCensoring was',dataOptionsStructure$censoringType,'with level sigma =',dataOptionsStructure$censoringLevel,fill=TRUE)
		if(dataOptionsStructure$covFuncFormGen=='Matern'){cat('\tCovariance function was',dataOptionsStructure$covFuncFormGen,'with Matern parameter set to',dataOptionsStructure$extraParamGen,fill=TRUE)
		} else {cat('\tCovariance function was',dataOptionsStructure$covFuncFormGen,fill=TRUE)}
		cat('\tMean function was',dataOptionsStructure$meanFuncFormGen,fill=TRUE)

		cat('Fitting options:',fill=TRUE)
		cat('\tPrediction was using aft',fill=TRUE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#---------------------------------------- Extract synthetic data ---------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	dimension 			= trainingTestStructure$dimension
	trainingAll 		= data.frame(trainingTestStructure$trainingData)
	names(trainingAll) 	= paste0('x',1:trainingTestStructure$dimension)
	trainingAll$y 		= trainingTestStructure$trainingTargets
	trainingAll$event 	= trainingTestStructure$events

	testAll 			= data.frame(trainingTestStructure$testData)
	names(testAll) 		= paste0('x',1:trainingTestStructure$dimension)
	if(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse'){
		testAll$y 		= trainingTestStructure$testTargets
		testAll$event 	= trainingTestStructure$testEvents
	} else if(dataOptionsStructure$censoringType=='None'){ 
		testAll$y 		= trainingTestStructure$testTargets
		testAll$event 	= rep(1,dim(trainingTestStructure$testTargets)[1]) 
	} else {
		testAll$y 		= trainingTestStructure$testTargetsPreCensoring
		testAll$event 	= rep(1,dim(trainingTestStructure$testTargets)[1])
	}

	timeStart 			= Sys.time()

	#-------------------------------------------------------------------------------------------------------#
	#----------------------------------- Accelerated Failure Time Model ------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	model.aft 			<- survreg(as.formula(paste0('Surv(y,event) ~',paste(paste0('x',1:dimension),collapse='+'))),data=trainingAll, dist="weibull")
	trainingPredictions <- predict(model.aft,type='response')
	testPredictions 	<- predict(model.aft,testAll,type='response')
	metricStructure  	<- CalculateMetrics(testPredictions,testAll$y,testAll$event,trainingPredictions,trainingAll$y,trainingAll$event)
	c.index 			<- metricStructure$c.index
	rmse 				<- metricStructure$rmse

	timeEnd   			<- Sys.time()
	timeTaken 			<- difftime(timeEnd,timeStart, units='min')

	#-------------------------------------------------------------------------------------------------------#
	#--------------------------------- Print and save variables and plots ----------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	if(printResults){
		cat('-------------------------------------------------------',fill=TRUE)
		cat('Model = ',model,fill=TRUE)
		cat('Time taken = ',timeTaken,fill=TRUE)
		cat('C Index = ',c.index,fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}

	if(printPlots){
		if(!any(is.infinite(testPredictions))){
			plot(testAll$y,testPredictions,xlim=c(0,max(c(testAll$y,testPredictions))),ylim=c(0,max(c(testAll$y,testPredictions))),ylab='Predicted Survival, AFT',xlab='Measured Survival',main=paste0('c.index = ',round(c.index,2)))
			abline(0,1)
		}
	}

	if(printPlots|savePlots){
		PlotKaplanMeier(testPredictions,testAll$y,testAll$event,model)
		plotKM <- recordPlot()
		# plotKM <- NULL
	}

	if(saveFiles){
		sink(paste0(fileName,'AFTMetrics.txt'),append=TRUE)
			cat('Model =',model,fill=TRUE)
			cat('Time taken = ',timeTaken,fill=TRUE)
			cat('C Index = ',c.index,fill=TRUE)
			cat('-------------------------------------------------------',fill=TRUE)
		sink()
		write.table(testPredictions,paste0(fileName,'TestTargetPredictions.csv'),sep=',',quote=FALSE,row.names=FALSE,append=TRUE)
		write.table(testAll$y,paste0(fileName,'testTargetsPreCensoring.csv'),sep=',',quote=FALSE,row.names=FALSE,append=TRUE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------------------- Return output --------------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	toReturn = list('timeTaken'=timeTaken,'trainingPredictions'=trainingPredictions,'testPredictions'=testPredictions,'c.index'=c.index,'rmse'=rmse,
					'model'=model,'trainingTestStructure'=trainingTestStructure)
	if(printPlots) toReturn$'plotKM' <- plotKM

	return(toReturn)
}
	