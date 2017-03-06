SetParametersRealYuan <- function(cancer,molPlatform,clinicalFlag,nReps,unid){
	#---------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#---------------------------------------------------------------------------------------#
	# Setting parameters for runRealYuanEtAl.R
	#---------------------------------------------------------------------------------------#

	# unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
	outerFolder 			<- 'Runs'
	folderName 				<- paste0(outerFolder,'/',unid)
	
	modelsList 				<- c('GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV')


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
									'printResults'=printResults,'folderName'=folderName,'outerFolder'=outerFolder)


	##-------------------------------------------------------------------------------------##
	##---------------------------------- Data Parameters ----------------------------------##
	##-------------------------------------------------------------------------------------##
	dataSource 				<- 'TCGAYuanEtAl' 		# 'R', 'Generate', 'RCensored', 'CuratedOvarian' (for CuratedOvarian attention needs to be paid to options below)
	dimension 				<- 10 				# 1 informative + 0 un-informative
	proportionTest 			<- NA
	nTraining 				<- NA
	nTest 					<- NA
	if(dataSource=='TCGAYuanEtAl'){
		nReps 				<- nReps
		cancer 				<- cancer
		molPlatform 		<- molPlatform
		clinicalFlag		<- clinicalFlag
	}
	
	logHypGenerate 			<- list('noise'=log(1.1),'func'=log(1.1),'length'=log(1.1),'mean'=c(rep(0,dimension),0))
	
	covFuncFormGen 			<- 'SqExp' 						# 'SqExp', 'Matern', 'ARD'
	extraParamGen 			<- 3
	meanFuncFormGen 		<- 'Linear'
	gridMinimum 			<- 0
	gridMaximum 			<- 8
	censoringType 			<- 'NormalLoopSample' 			# 'Normal', 'None' or 'NormalLoopSample' ('NormalLoopRerun' is not non-informative censoring)
	censoringSD 			<- 50
	censoringMean 			<- 0
	censoredProportion 		<- 50
	nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)
	if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE
	extraDimensions 		<- 0 

	dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,'meanFuncFormGen'=meanFuncFormGen,
									'extraParamGen'=extraParamGen,'dimension'=dimension,'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,
									'gridMaximum'=gridMaximum,'censoringSD'=censoringSD,'censoringMean'=censoringMean,'censoringType'=censoringType,'nCensored'=nCensored,
									'useARD'=useARD,'extraDimensions'=extraDimensions,'folderName'=folderName,'proportionTest'=proportionTest,
									'censoredProportion'=censoredProportion,'nReps'=nReps,'cancer'=cancer,'molPlatform'=molPlatform,'clinicalFlag'=clinicalFlag)

	##-------------------------------------------------------------------------------------##
	##---------------------------- Gaussian Process Parameters ----------------------------##
	##-------------------------------------------------------------------------------------##
	modelType 				<- 'survival' 			# 'non-survival' or 'survival'
	covFuncForm 			<- 'SqExp'				# 'SqExp', 'Matern', 'ARD', 'InformedARD', 'RF'
	meanFuncForm 			<- 'Zero'
	burnIn 					<- FALSE
	extraParam 				<- 3
	tolerance 				<- 0.0001*nTraining*censoredProportion
	hypChangeTol 			<- 10^-10
	maxCount 				<- 500
	maxit 					<- 2000
	maxitPreLearn 			<- 100
	maxitSurvival 			<- 2
	maxitLaplace			<- 2
	optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
	noiseCorr 				<- FALSE
	imposePriors 			<- TRUE

	logHypStart 			<- list('noise'=log(0.01),'func'=log(1.1),'length'=log(0.5),'mean'=c(rep(0,dimension),0))

	parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'extraParam'=extraParam,'maxit'=maxit,'maxitPreLearn'=maxitPreLearn,
									'maxitSurvival'=maxitSurvival,'optimType'=optimType,'logHypStart'=logHypStart,'modelType'=modelType,'unid'=unid,'tolerance'=tolerance,
									'maxCount'=maxCount,'hypChangeTol'=hypChangeTol,'burnIn'=burnIn,'maxitLaplace'=maxitLaplace,'noiseCorr'=noiseCorr,
									'imposePriors'=imposePriors)

	if(saveFiles){
		dir.create(file.path(paste0(getwd(),"/",outerFolder),unid),showWarnings=FALSE)
		if(dataSource=='Generate'){
			dir.create(file.path(paste0(getwd(),"/",outerFolder,'/',unid),'data'),showWarnings=FALSE)
			write.table(exp(unlist(logHypStart))[1:3],paste0(folderName,'/data/hypStart.csv'),row.names=FALSE,col.names=FALSE,sep=',')
			write.table(exp(unlist(logHypGenerate))[1:3],paste0(folderName,'/data/hypGenerate.csv'),row.names=FALSE,col.names=FALSE,sep=',')
		}
	}

	# if(saveFiles){
	# 	sink(paste0(getwd(),'/',outerFolder,'/',unid,'/',unid,'RunOptions.txt'),append=TRUE)
	# 		PrintOptions2(dataOptionsStructure,parameterStructure)
	# 	sink()
	# }

	toReturn <- list('plotSaveOptions'=plotSaveOptions,'dataOptionsStructure'=dataOptionsStructure,'parameterStructure'=parameterStructure)

}