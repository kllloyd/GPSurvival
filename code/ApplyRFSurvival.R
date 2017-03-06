ApplyRFSurvival <- function(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions){
    #-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Predicting risk of data set members given censored and uncensored training data.
	# Applied Random Forest model for survival data using rfsrc
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

	model 			= 'RFSurvival'
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
		cat('\tPrediction was using random forest for survival data, rfsrc',fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#---------------------------------------- Extract synthetic data ---------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	dimension 			= trainingTestStructure$dimension
	trainingAll 		= data.frame(trainingTestStructure$trainingData)
	names(trainingAll) 	= paste0('x',1:trainingTestStructure$dimension)
	trainingAll$V1 		= c(trainingTestStructure$trainingTargets)
	trainingAll$V1.1 	= c(trainingTestStructure$events)

	testAll 			= data.frame(trainingTestStructure$testData)
	names(testAll) 		= paste0('x',1:trainingTestStructure$dimension)
	if(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse'|dataOptionsStructure$dataSource=='TCGAYuanEtAl'){
		testAll$V1 		= c(trainingTestStructure$testTargets)
		testAll$V1.1 	= c(trainingTestStructure$testEvents)
	} else if(dataOptionsStructure$censoringType=='None'){ 
		testAll$V1 		= c(trainingTestStructure$testTargets)
		testAll$V1.1 	= c(rep(1,dim(trainingTestStructure$testTargets)[1]))
	} else {
		testAll$V1 		= c(trainingTestStructure$testTargetsPreCensoring)
		testAll$V1.1 	= c(rep(1,dim(trainingTestStructure$testTargets)[1]))
	}

	timeStart 			= Sys.time()
	
	#-------------------------------------------------------------------------------------------------------#
	#---------------------------------------- Random Forest Model ------------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	model.rfSurvival 	<- rfsrc(Surv(V1,V1.1)~.,data=trainingAll,na.action='na.impute',importance=TRUE)
	modelCoefficients 	<- model.rfSurvival$importance
	predTrain 			<- predict(model.rfSurvival)
	trainingPredictions <- predTrain$predicted
	predTest 			<- predict(model.rfSurvival,testAll)
	survivalFunc 		<- predTest$survival
	timeFunc 			<- predTest$time.interest
	testPredictions 	<- predTest$predicted
	metricStructure  	<- CalculateMetrics(testPredictions,testAll$V1,testAll$V1.1,trainingPredictions,trainingAll$V1,trainingAll$V1.1)
	c.index 			<- metricStructure$c.index
	rmse 				<- metricStructure$rmse

	timeEnd   			<- Sys.time()
	timeTaken 			<- difftime(timeEnd,timeStart, units='min')
	genesByImportance 	<- colnames(trainingTestStructure$trainingData)[order(model.rfSurvival$importance)]
	# browser()
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
		plot(testPredictions,testAll$V1,col='white',xlim=c(0,max(c(testAll$V1,testPredictions))),ylim=c(0,max(c(testAll$V1,testPredictions))),xlab='Predicted Survival, Random Forest for Survival',ylab='Measured Survival',main=paste0('c.index = ',round(c.index,2)))
		points(testPredictions[testAll$V1.1==1],testAll$V1[testAll$V1.1==1],col='blue')
		points(testPredictions[testAll$V1.1==0],testAll$V1[testAll$V1.1==0],col='green')
		abline(0,1)

		PlotKaplanMeier(testPredictions,testAll$V1,testAll$V1.1,model)
		plotKM <- recordPlot()
	} else {
		plotKM <- NULL
	}

	if(saveFiles){
		sink(paste0(fileName,'RFSurvivalMetrics.txt'),append=TRUE)
			cat('Model = RFSurvival',fill=TRUE)
			cat('Time taken = ',timeTaken,fill=TRUE)
			cat('C Index = ',c.index,fill=TRUE)
			cat('-------------------------------------------------------',fill=TRUE)
		sink()
		write.table(testPredictions,paste0(fileName,'RFSurvivalTestTargetPredictions.csv'),sep=',',quote=FALSE,row.names=FALSE)
		write.table(testAll$V1,paste0(fileName,'RFSurvivalTestTargets.csv'),sep=',',quote=FALSE,row.names=FALSE)
	}

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------------------- Return output --------------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	toReturn = list('timeTaken'=timeTaken,'trainingPredictions'=trainingPredictions,'testPredictions'=testPredictions,'c.index'=c.index,'rmse'=rmse,
					'model'='RFSurvival','modelVariables'=modelCoefficients,'genesByImportance'=genesByImportance,'testTargets'=trainingTestStructure$testTargets,
					'trainingTestStructure'=trainingTestStructure)
	if(printPlots) toReturn$'plotKM' <- plotKM

	return(toReturn)
}
	