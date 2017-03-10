SetParametersExp3 <- function(nTraining,nTest,hypGenerateNoise,censoredProportion,nReps,unid){
	#---------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#---------------------------------------------------------------------------------------#
	# Setting parameters for runSyntheticExp3.R
	#---------------------------------------------------------------------------------------#

	# unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
	outerFolder 			<- 'Runs'
	folderName 				<- paste0(outerFolder,'/',unid)
	
	# nReps 					<- nReps
	
	modelsList 				<- c('GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV','GPSurvLaplace','AFT','Coxph','Glmnet','GBM','RF','RFSurvival')


	##-------------------------------------------------------------------------------------##
	##----------------------------- Plot and Save Parameters ------------------------------##
	##-------------------------------------------------------------------------------------##
	printOptions 			<- FALSE
	printPlots 				<- TRUE
	savePlots 				<- TRUE
	saveFiles 				<- TRUE
	plotRuns 				<- TRUE
	printResults 			<- TRUE
	plotSaveOptions 		<- list('printOptions'=printOptions,'printPlots'=printPlots,'savePlots'=savePlots,'saveFiles'=saveFiles,'plotRuns'=plotRuns,
									'printResults'=printResults,'folderName'=folderName,'outerFolder'=outerFolder,'nReps'=nReps)


	##-------------------------------------------------------------------------------------##
	##---------------------------------- Data Parameters ----------------------------------##
	##-------------------------------------------------------------------------------------##
	dataSource 				<- 'Generate' 		# 'R', 'Generate', 'RCensored', 'CuratedOvarian' (for CuratedOvarian attention needs to be paid to options below)
	dimension 				<- 6 				# 1 informative + 0 un-informative
	proportionTest 			<- NA
	# nTraining 				<- nTraining
	# nTest 					<- nTest 
	
	logHypGenerate 			<- list('noise'=log(hypGenerateNoise),'func'=log(0.5),'length'=log(1.1),'mean'=c(rep(0,dimension),0))
	
	covFuncFormGen 			<- 'SqExp'
	extraParamGen 			<- 3
	meanFuncFormGen 		<- 'Linear'
	gridMinimum 			<- 0
	gridMaximum 			<- 8
	censoringType 			<- 'NormalLoopSample' 			# 'Normal', 'None' or 'NormalLoopSample' ('NormalLoopRerun' is not non-informative censoring)
	censoringSD 			<- 50
	censoringMean 			<- 0
	censoredProportion 		<- censoredProportion
	nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)
	if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE
	extraDimensions 		<- 0 

	if(dataSource=='CuratedOvarian'){
		geneSubsetFlag 			<- 'SystematicReview' 										# 'All', 'SystematicReview', 'TaqmanArray', 'MP3'
		geneExpressionFlag 		<- TRUE
		clinicalFeaturesFlag 	<- FALSE
		dataSet 				<- 'GSE9891'
	}

	dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,'meanFuncFormGen'=meanFuncFormGen,
									'extraParamGen'=extraParamGen,'dimension'=dimension,'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,
									'gridMaximum'=gridMaximum,'censoringSD'=censoringSD,'censoringMean'=censoringMean,'censoringType'=censoringType,'nCensored'=nCensored,
									'useARD'=useARD,'extraDimensions'=extraDimensions,'folderName'=folderName,'proportionTest'=proportionTest,
									'censoredProportion'=censoredProportion)

	if(dataSource=='CuratedOvarian'){
		dataOptionsStructure$chosenClinicalFeatures 		<- chosenClinicalFeatures
		dataOptionsStructure$chosenExtraClinicalFeatures 	<- chosenExtraClinicalFeatures
		dataOptionsStructure$chosenExpressionFeatures 		<- chosenExpressionFeatures
		dataOptionsStructure$geneSubsetFlag 				<- geneSubsetFlag
		dataOptionsStructure$geneExpressionFlag 			<- geneExpressionFlag
		dataOptionsStructure$clinicalFeaturesFlag 			<- clinicalFeaturesFlag
		dataOptionsStructure$dataSet 						<- dataSet
	}


	##-------------------------------------------------------------------------------------##
	##---------------------------- Gaussian Process Parameters ----------------------------##
	##-------------------------------------------------------------------------------------##
	modelType 				<- 'survival' 			# 'non-survival' or 'survival'
	covFuncForm 			<- 'SqExp'				# 'SqExp' or 'ARD'
	meanFuncForm 			<- 'Zero'
	burnIn 					<- FALSE
	extraParam 				<- 3
	tolerance 				<- 0.00001*nTraining*censoredProportion
	hypChangeTol 			<- 10^-10
	maxCount 				<- 500
	maxit 					<- 2000
	maxitPreLearn 			<- 1
	maxitSurvival 			<- 2
	maxitLaplace			<- 2
	optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
	noiseCorr 				<- FALSE
	imposePriors 			<- TRUE
	if(dataSource=='Generate'){
		logHypStart 	<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(1.2),'mean'=c(rep(0,dimension),0))
	} else if (datatSource=='CuratedOvarian'){
		forHyp 			= matrix(log(abs(runif(1+1+dimension+dimension+1,min=0.005,max=5))),nrow=1)		# starting hyperparameters drawn from uniform distribution
		if(covFuncForm=='ARD'){
			logHypStart <- list('noise'=forHyp[1],'func'=forHyp[2],'length'=forHyp[3:(dimension+2)],'mean'=exp(forHyp[(dim(forHyp)[2]-dimension):dim(forHyp)[2]]))
		} else {
			logHypStart <- list('noise'=forHyp[1],'func'=forHyp[2],'length'=forHyp[3],'mean'=exp(forHyp[(dim(forHyp)[2]-dimension):dim(forHyp)[2]]))
		}
	} 
	parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'extraParam'=extraParam,'maxit'=maxit,'maxitPreLearn'=maxitPreLearn,
									'maxitSurvival'=maxitSurvival,'optimType'=optimType,'logHypStart'=logHypStart,'modelType'=modelType,'unid'=unid,'tolerance'=tolerance,
									'maxCount'=maxCount,'hypChangeTol'=hypChangeTol,'burnIn'=burnIn,'maxitLaplace'=maxitLaplace,'noiseCorr'=noiseCorr,
									'imposePriors'=imposePriors)

	if(printOptions) PrintOptions(dataOptionsStructure,parameterStructure,modelsList)
	if(saveFiles){
		dir.create(file.path(paste0(getwd(),"/",outerFolder),unid),showWarnings=FALSE)
		if(dataSource=='Generate'){
			dir.create(file.path(paste0(getwd(),"/",outerFolder,'/',unid),'data'),showWarnings=FALSE)
			write.table(exp(unlist(logHypStart))[1:3],paste0(folderName,'/data/hypStart.csv'),row.names=FALSE,col.names=FALSE,sep=',')
			write.table(exp(unlist(logHypGenerate))[1:3],paste0(folderName,'/data/hypGenerate.csv'),row.names=FALSE,col.names=FALSE,sep=',')
		}
		sink(paste0(getwd(),'/',outerFolder,'/',unid,'/',unid,'RunOptions.txt'))
		PrintOptions2(dataOptionsStructure,parameterStructure)
		sink()
	}

	toReturn <- list('plotSaveOptions'=plotSaveOptions,'dataOptionsStructure'=dataOptionsStructure,'parameterStructure'=parameterStructure)

}