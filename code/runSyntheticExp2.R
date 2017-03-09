#----------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
# Applying Gaussian process to synthetic data generated with censored survival times in 
# the training set.
# All models applied to the same data. Models are:	GP, GPS1, GPS2, GPS3, AFT, Cox PH, 
# 													Glmnet, GBM and RSF
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Boxplot to compare concordance index results of different models
#---------------------------------------------------------------------------------------#

##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##

library(fields)
library(gbm)
library(glmnet)
library(ipred)
library(MASS)
library(Matrix)
library(nlme)
library(NORMT3)
library(pdist)
library(randomForestSRC)
library(rgl)
library(rms)
library(survcomp)
library(survival)
library(zoo)


##-------------------------------------------------------------------------------------##
##--------------------------------- Source Functions ----------------------------------##
##-------------------------------------------------------------------------------------##
source('add.alpha.R')
source('AdjustTrainingSurvivalMeanVariance.R')
source('ApplyAFT.R')
source('ApplyCoxph.R')
source('ApplyGBM.R')
source('ApplyGlmnet.R')
source('ApplyGP.R')
source('ApplyGPDiffInf.R')
source('ApplyRF.R')
source('ApplyRFSurvival.R')
source('CalculateMetrics.R')
source('CensorData.R')
source('CovFunc.R')
source('DataExtraction.R')
source('GenerateData.R')
source('ImputeMissingData.R')
source('LogPriorX.R')
source('MakeSyntheticData.R')
source('MeanFunc.R')
source('NormaliseExpressionData.R')
source('PlotKaplanMeier.R')
source('PreLearnHyperparam.R')
source('PrintOptions.R')
source('PrintOptions2.R')
source('RemoveCensored.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##

set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
outerFolder 			<- 'Runs'
folderName 				<- paste0(outerFolder,'/',unid)
dir.create(file.path(getwd(),outerFolder),showWarnings=FALSE)

nReps 					<- 30

modelsList 				<- c('GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV','GPSurvInfUnif','GPSurvInfMed','AFT','Coxph','Glmnet','GBM','RF','RFSurvival')


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
dataSource 				<- 'Generate' 
dimension 				<- 4 
extraDimensions 		<- 1 
proportionTest 			<- NA
nTraining 				<- 200
nTest 					<- 50 

covFuncFormGen 			<- 'ARD'
maternParamGen 			<- 3
meanFuncFormGen 		<- 'Linear'
if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE

if(useARD){
	logHypGenerate 		<- list('noise'=log(0.01),'func'=log(0.5),'length'=log(c(1.1,0.9,1.2,0.7)),'mean'=c(rep(0,dimension),0))
} else {
	logHypGenerate 		<- list('noise'=log(0.01),'func'=log(0.5),'length'=log(1.1),'mean'=c(rep(0,dimension),0))
}

gridMinimum 			<- 0
gridMaximum 			<- 8
censoringType 			<- 'NormalLoopSample'
censoringSD 			<- 50
censoringMean 			<- 0
censoredProportion 		<- 0.70
nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)

dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,'meanFuncFormGen'=meanFuncFormGen,
								'maternParamGen'=maternParamGen,'dimension'=dimension,'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,
								'gridMaximum'=gridMaximum,'censoringSD'=censoringSD,'censoringMean'=censoringMean,'censoringType'=censoringType,'nCensored'=nCensored,
								'useARD'=useARD,'extraDimensions'=extraDimensions,'folderName'=folderName,'proportionTest'=proportionTest,
								'censoredProportion'=censoredProportion)


##-------------------------------------------------------------------------------------##
##---------------------------- Gaussian Process Parameters ----------------------------##
##-------------------------------------------------------------------------------------##
modelType 				<- 'survival' 			# 'non-survival' or 'survival'
covFuncForm 			<- 'SqExp'
meanFuncForm 			<- 'Zero'
burnIn 					<- FALSE
maternParam 			<- 3
tolerance 				<- 0.0001*nTraining*censoredProportion
toleranceLaplace 		<- 0.0001*nTraining*censoredProportion
hypChangeTol 			<- 10^-10
maxCount 				<- 500
maxit 					<- 2000
maxitPreLearn 			<- 1
maxitSurvival 			<- 10
maxitLaplaceHyp			<- 2
maxitLaplaceFHat 		<- 100
optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
noiseCorr 				<- FALSE
imposePriors 			<- FALSE
if(!useARD){
	logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dimension),0))
} else {
	logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=rep(log(0.9),dimension),'mean'=c(rep(0,dimension),0))

}
parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'maternParam'=maternParam,'maxit'=maxit,'maxitPreLearn'=maxitPreLearn,
								'maxitSurvival'=maxitSurvival,'optimType'=optimType,'logHypStart'=logHypStart,'modelType'=modelType,'unid'=unid,'tolerance'=tolerance,
								'toleranceLaplace'=toleranceLaplace,'maxCount'=maxCount,'hypChangeTol'=hypChangeTol,'burnIn'=burnIn,'maxitLaplaceHyp'=maxitLaplaceHyp,
								'maxitLaplaceFHat'=maxitLaplaceFHat,'noiseCorr'=noiseCorr,'imposePriors'=imposePriors)

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


##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureGPNonSurvNoCens 	<- list()
outputStructureGPSurvNoCorr 	<- list()
outputStructureGPSurvCorrV		<- list()
outputStructureGPSurvCorrL 		<- list()
outputStructureGPSurvInfUnif 	<- list()
outputStructureGPSurvInfMed 	<- list()
outputStructureAFT				<- list()
outputStructureCoxph 			<- list()
outputStructureGBM				<- list()
outputStructureGlmnet 			<- list()
outputStructureRF				<- list()
outputStructureRFSurvival 		<- list()

c.index.GPNonSurvNoCens 		<- rep(NA,nReps)
c.index.GPSurvNoCorr 			<- rep(NA,nReps)
c.index.GPSurvCorrV				<- rep(NA,nReps)
c.index.GPSurvCorrL 			<- rep(NA,nReps)
c.index.GPSurvInfUnif 			<- rep(NA,nReps)
c.index.GPSurvInfMed 			<- rep(NA,nReps)
c.index.AFT						<- rep(NA,nReps)
c.index.Coxph 					<- rep(NA,nReps)
c.index.GBM						<- rep(NA,nReps)
c.index.Glmnet 					<- rep(NA,nReps)
c.index.RF						<- rep(NA,nReps)
c.index.RFSurvival 				<- rep(NA,nReps)

rmse.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.GPSurvNoCorr 				<- rep(NA,nReps)
rmse.GPSurvCorrV				<- rep(NA,nReps)
rmse.GPSurvCorrL 				<- rep(NA,nReps)
rmse.GPSurvInfUnif 				<- rep(NA,nReps)
rmse.GPSurvInfMed 				<- rep(NA,nReps)
rmse.AFT						<- rep(NA,nReps)
rmse.Coxph 						<- rep(NA,nReps)
rmse.GBM						<- rep(NA,nReps)
rmse.Glmnet 					<- rep(NA,nReps)
rmse.RF							<- rep(NA,nReps)
rmse.RFSurvival 				<- rep(NA,nReps)


##-------------------------------------------------------------------------------------##
##----------------------------------- Generate Data -----------------------------------##
##-------------------------------------------------------------------------------------##
trainingTestStructureForNReps 			<- GenerateData(dataOptionsStructure,outerFolder,nReps)
for(i in 1:nReps){
	trainingTestStructureForNReps[[i]] 	<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
	cat('Data set for run',i,'normalised',fill=TRUE)
}


##-------------------------------------------------------------------------------------##
##----------------------------- Run GPNonSurvNoCens Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureGPNonSurvNoCens[[i]] <- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$c.index)!=0,outputStructureGPNonSurvNoCens[[i]]$c.index,NA)
	rmse.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$rmse)!=0,outputStructureGPNonSurvNoCens[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS1 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvNoCorr[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvNoCorr[i] 			<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$c.index)!=0,outputStructureGPSurvNoCorr[[i]]$c.index,NA)
	rmse.GPSurvNoCorr[i] 				<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$rmse)!=0,outputStructureGPSurvNoCorr[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS2 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'noiseCorrLearned'
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvCorrL[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$c.index)!=0,outputStructureGPSurvCorrL[[i]]$c.index,NA)
	rmse.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$rmse)!=0,outputStructureGPSurvCorrL[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS3 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvCorrV[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$c.index)!=0,outputStructureGPSurvCorrV[[i]]$c.index,NA)
	rmse.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$rmse)!=0,outputStructureGPSurvCorrV[[i]]$rmse,NA)
}

##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPR1 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'none'
parameterStructure$modelType 			<- 'survival'
parameterStructure$inferenceType 		<- 'median'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInfMed[[i]] 	<- ApplyGPDiffInf(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInfMed[i] 			<- ifelse(length(outputStructureGPSurvInfMed[[i]]$c.index)!=0,outputStructureGPSurvInfMed[[i]]$c.index,NA)
	rmse.GPSurvInfMed[i] 				<- ifelse(length(outputStructureGPSurvInfMed[[i]]$rmse)!=0,outputStructureGPSurvInfMed[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPR2 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'none'
parameterStructure$modelType 			<- 'survival'
parameterStructure$inferenceType 		<- 'uniform'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInfUnif[[i]] 	<- ApplyGPDiffInf(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInfUnif[i] 			<- ifelse(length(outputStructureGPSurvInfUnif[[i]]$c.index)!=0,outputStructureGPSurvInfUnif[[i]]$c.index,NA)
	rmse.GPSurvInfUnif[i] 				<- ifelse(length(outputStructureGPSurvInfUnif[[i]]$rmse)!=0,outputStructureGPSurvInfUnif[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##----------------------------------- Run AFT Model -----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureAFT[[i]] 			<- ApplyAFT(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.AFT[i] 						<- ifelse(length(outputStructureAFT[[i]]$c.index)!=0,outputStructureAFT[[i]]$c.index,NA)
	rmse.AFT[i] 						<- ifelse(length(outputStructureAFT[[i]]$rmse)!=0,outputStructureAFT[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Coxph Model ----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureCoxph[[i]] 			<- ApplyCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Coxph[i] 					<- ifelse(length(outputStructureCoxph[[i]]$c.index)!=0,outputStructureCoxph[[i]]$c.index,NA)
	rmse.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$rmse)!=0,outputStructureCoxph[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##----------------------------------- Run GBM Model -----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGBM[[i]] 			<- ApplyGBM(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GBM[i] 						<- ifelse(length(outputStructureGBM[[i]]$c.index)!=0,outputStructureGBM[[i]]$c.index,NA)
	rmse.GBM[i] 						<- ifelse(length(outputStructureGBM[[i]]$rmse)!=0,outputStructureGBM[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Glmnet Model ---------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGlmnet[[i]] 			<- ApplyGlmnet(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Glmnet[i] 					<- ifelse(length(outputStructureGlmnet[[i]]$c.index)!=0,outputStructureGlmnet[[i]]$c.index,NA)
	rmse.Glmnet[i] 						<- ifelse(length(outputStructureGlmnet[[i]]$rmse)!=0,outputStructureGlmnet[[i]]$rmse,NA)
}


#-------------------------------------------------------------------------------------##
#------------------------------------ Run RF Model -----------------------------------##
#-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureRF[[i]] 				<- ApplyRF(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.RF[i] 						<- ifelse(length(outputStructureRF[[i]]$c.index)!=0,outputStructureRF[[i]]$c.index,NA)
	rmse.RF[i] 							<- ifelse(length(outputStructureRF[[i]]$rmse)!=0,outputStructureRF[[i]]$rmse,NA)
}


#-------------------------------------------------------------------------------------##
#------------------------------- Run RFSurvival Model --------------------------------##
#-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureRFSurvival[[i]] 		<- ApplyRFSurvival(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.RFSurvival[i] 				<- ifelse(length(outputStructureRFSurvival[[i]]$c.index)!=0,outputStructureRFSurvival[[i]]$c.index,NA)
	rmse.RFSurvival[i] 					<- ifelse(length(outputStructureRFSurvival[[i]]$rmse)!=0,outputStructureRFSurvival[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------ Print and Save Results -------------------------------##
##-------------------------------------------------------------------------------------##
							##-------------------------------##
							##---- Print C Index & RMSE -----##
							##-------------------------------##
if(printResults){
	cat('---------------------------------------',fill=TRUE)
	cat('Concordance Indices:',fill=TRUE)
	cat('GPNonSurvNoCens mean c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS1 mean c index =',paste0(round(c.index.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPS2 mean c index =',paste0(round(c.index.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('GPS3 mean c index =',paste0(round(c.index.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('AFT mean c index =',paste0(round(c.index.AFT,4),collapse=', '),fill=TRUE)
	cat('Coxph mean c index =',paste0(round(c.index.Coxph,4),collapse=', '),fill=TRUE)
	cat('GBM mean c index =',paste0(round(c.index.GBM,4),collapse=', '),fill=TRUE)
	cat('Glmnet mean c index =',paste0(round(c.index.Glmnet,4),collapse=', '),fill=TRUE)
	cat('RF mean c index =',paste0(round(c.index.RF,4),collapse=', '),fill=TRUE)
	cat('RFSurvival mean c index =',paste0(round(c.index.RFSurvival,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
	cat('RMSE:',fill=TRUE)
	cat('GP mean rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS1 mean rmse =',paste0(round(rmse.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPS2 mean rmse =',paste0(round(rmse.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('GPS3 mean rmse =',paste0(round(rmse.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
}

							##-------------------------------##
							##--- Plot Kaplan-Meier Plots ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvNoCorr[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrL[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrV[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvInfUnif[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvInfMed[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','AFT','/',unid,'AFT','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureAFT[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Coxph','/',unid,'Coxph','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureCoxph[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GBM','/',unid,'GBM','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGBM[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Glmnet','/',unid,'Glmnet','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGlmnet[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RF','/',unid,'RF','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureRF[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RFSurvival','/',unid,'RFSurvival','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureRFSurvival[[i]]$plotKM)
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot Measured/Predicted ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPNonSurvNoCens[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvNoCorr[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrL[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrV[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvInfUnif[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvInfMed[[i]]$plot2)
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot GP Hyperparameters ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(rbind(unlist(t(outputStructureGPNonSurvNoCens[[i]]$logHypChosen[1:3])),unlist(t(outputStructureGPNonSurvNoCens[[i]]$parameterStructure$logHypStart[1:3])))
					/unlist(outputStructureGPNonSurvNoCens[[i]]$dataOptionsStructure$logHypGenerate)[1:3],
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvNoCorr[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvNoCorr[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvCorrL[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvCorrL[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvCorrV[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvCorrV[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvInfUnif[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvInfUnif[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvInfMed[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvInfMed[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)
}


							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##
modelNames 		<- c('AFT','Coxph','Glmnet','GBM','RFSurvival','GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV','GPSurvInfMed','GPSurvInfUnif') # Reordered version of modelsList
modelNamesPlot 	<- c('AFT','Cox PH','Glmnet','GBM','RSF','GP','GPS1','GPS2','GPS3','GPR1','GPR2')
toInvert 		<- c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
for(i in 1:length(modelNames)){
	c.index.mat[,i] 	<- get(paste0('c.index.',modelNames[i]))
	if(toInvert[i]) c.index.mat[,i] <- 1-c.index.mat[,i]
	c.index.mean[i] 	<- mean(c.index.mat[,i],na.rm=TRUE)
}

if(savePlots){
	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotCIndexAllModels.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		plot(sort(rep(1:length(modelNames),nReps)),c(c.index.mat),pch=20,col=add.alpha('chartreuse3',0.8),cex=0.7,xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1))
		points(1:length(modelNames),c.index.mean,pch=20,cex=1.2)
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		text(1:length(modelNames),rep(0.97,length(modelNames)),labels=round(c.index.mean,4),cex=0.7,pos=3)
		layout(1)
	dev.off(pdf.output)

	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotCIndexAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(c.index.mat,col=add.alpha('chartreuse3',0.8),xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1))
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}


##-------------------------------------------------------------------------------------##
##------------------------------------ Save Output ------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureAll <- list('outputStructureGPNonSurvNoCens'=outputStructureGPNonSurvNoCens,
							'outputStructureGPSurvCorrV'=outputStructureGPSurvCorrV,
							'outputStructureRF'=outputStructureRF,
							'outputStructureRFSurvival'=outputStructureRFSurvival,
							'outputStructureGPSurvCorrL'=outputStructureGPSurvCorrL,
							'outputStructureAFT'=outputStructureAFT,
							'outputStructureCoxph'=outputStructureCoxph,
							'outputStructureGlmnet'=outputStructureGlmnet,
							'outputStructureGBM'=outputStructureGBM)

save(list='outputStructureAll',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureAll','_','Workspace.RData'))