SetParametersRealTothill <- function(geneSubsetFlag,clinicalSubsetFlag,dimension,nReps){
	#---------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#---------------------------------------------------------------------------------------#
	# Setting parameters for runRealTothillEtAl.R
	#---------------------------------------------------------------------------------------#

	unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
	outerFolder 			<- 'Runs'
	folderName 				<- paste0(outerFolder,'/',unid)

	nReps 					<- nReps

	modelsList 				<- c('AFT','Coxph','Glmnet','RFSurvival','GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV')


	##-------------------------------------------------------------------------------------##
	##----------------------------- Plot and Save Parameters ------------------------------##
	##-------------------------------------------------------------------------------------##
	printOptions 			<- FALSE
	printPlots 				<- TRUE
	savePlots 				<- TRUE
	saveFiles 				<- TRUE
	plotRuns 				<- TRUE
	printResults 			<- TRUE
	plotSaveOptions 		<- list('printOptions'=printOptions,'printPlots'=printPlots,'savePlots'=savePlots,'saveFiles'=saveFiles,
									'plotRuns'=plotRuns,'printResults'=printResults,'folderName'=folderName,'outerFolder'=outerFolder,
									'nReps'=nReps)


	##-------------------------------------------------------------------------------------##
	##---------------------------------- Data Parameters ----------------------------------##
	##-------------------------------------------------------------------------------------##
	dataSource 				<- 'CuratedOvarian'
	dimension 				<- dimension 		# GSE9891: clinical = 3, OCGS = 97, OCGS + clinical = 100, SRGS = 84, SRGS + clinical = 87
	proportionTest 			<- NA
	nTraining 				<- 200
	nTest 					<- 50

	logHypGenerate 			<- list('noise'=log(0.01),'func'=log(1.3),'length'=log(2.1),'mean'=c(rep(0,dimension),0))

	covFuncFormGen 			<- 'SqExp'
	maternParamGen 			<- 3
	meanFuncFormGen 		<- 'Linear'
	gridMinimum 			<- 0
	gridMaximum 			<- 8
	censoringType 			<- 'NormalLoopSample'
	censoringSD 			<- 50
	censoringMean 			<- 950
	censoredProportion 		<- 0.5
	nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)
	if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE
	extraDimensions 		<- 0 

	if(dataSource=='CuratedOvarian'){
		dataSet 					<- 'GSE9891'
		geneSubsetFlag 				<- geneSubsetFlag 								# 'All' = all genes, 'SRGS1' = SRGS, 'TaqMan' = OVCA, 'None' = no genes
		geneExpressionFlag 			<- ifelse(geneSubsetFlag=='None',FALSE,TRUE)
		clinicalSubsetFlag 			<- clinicalSubsetFlag							# 'One', 'Two', 'Three', etc. See GetDataSetFeatures.R
		clinicalFeaturesFlag 		<- c(NA,NA)
		clinicalFeaturesFlag[1] 	<- any(ifelse(clinicalSubsetFlag==c('One','Three','Four','Five','Six'),TRUE,FALSE))
		clinicalFeaturesFlag[2] 	<- any(ifelse(clinicalSubsetFlag==c('Two','Three','Six'),TRUE,FALSE))
		dataSetFeaturesStructure 	<- GetDataSetFeatures(dataSource,dataSet,geneSubsetFlag,clinicalSubsetFlag)
	}

	dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,
									'meanFuncFormGen'=meanFuncFormGen,'maternParamGen'=maternParamGen,'dimension'=dimension,
									'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,'gridMaximum'=gridMaximum,
									'censoringSD'=censoringSD,'censoringMean'=censoringMean,'censoringType'=censoringType,
									'nCensored'=nCensored,'useARD'=useARD,'extraDimensions'=extraDimensions,'folderName'=folderName,
									'proportionTest'=proportionTest,'censoredProportion'=censoredProportion)

	if(dataSource=='CuratedOvarian'){
		dataOptionsStructure$chosenClinicalFeatures 		<- dataSetFeaturesStructure$chosenClinicalFeatures
		dataOptionsStructure$chosenExpressionFeatures 		<- dataSetFeaturesStructure$chosenExpressionFeatures
		dataOptionsStructure$geneSubsetFlag 				<- geneSubsetFlag
		dataOptionsStructure$geneExpressionFlag 			<- geneExpressionFlag
		dataOptionsStructure$clinicalFeaturesFlag 			<- clinicalFeaturesFlag
		dataOptionsStructure$dataSet 						<- dataSet
		dataOptionsStructure$chosenExtraClinicalFeatures 	<- dataSetFeaturesStructure$chosenExtraClinicalFeatures
	}


	##-------------------------------------------------------------------------------------##
	##---------------------------- Gaussian Process Parameters ----------------------------##
	##-------------------------------------------------------------------------------------##
	modelType 				<- 'survival' 			# 'non-survival' or 'survival'
	covFuncForm 			<- 'SqExp'				# 'SqExp' or 'ARD'
	meanFuncForm 			<- 'Zero'
	burnIn 					<- FALSE
	maternParam 			<- 3
	tolerance 				<- 0.00001*nTraining*censoredProportion
	toleranceLaplace 		<- 0.00001*nTraining*censoredProportion
	hypChangeTol 			<- 10^-10
	maxCount 				<- 500
	maxit 					<- 2000
	maxitPreLearn 			<- 1
	maxitSurvival 			<- 2
	maxitLaplaceHyp			<- 2
	maxitLaplaceFHat 		<- 100
	optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
	noiseCorr 				<- FALSE
	imposePriors 			<- TRUE
	logHypStart 			<- list('noise'=log(0.9),'func'=log(0.9),'length'=log(0.9),'mean'=c(rep(0,dimension),0))

	parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'maternParam'=maternParam,'maxit'=maxit,
									'maxitPreLearn'=maxitPreLearn,'maxitSurvival'=maxitSurvival,'optimType'=optimType,
									'logHypStart'=logHypStart,'modelType'=modelType,'unid'=unid,'tolerance'=tolerance,
									'toleranceLaplace'=toleranceLaplace,'maxCount'=maxCount,'hypChangeTol'=hypChangeTol,'burnIn'=burnIn,
									'maxitLaplaceHyp'=maxitLaplaceHyp,'maxitLaplaceFHat'=maxitLaplaceFHat,'noiseCorr'=noiseCorr,
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