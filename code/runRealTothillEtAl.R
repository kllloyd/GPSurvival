#--------------------------------------------------------------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
#--------------------------------------------------------------------------------------------------------------------------------------------#
# Applying models to data extracted from package CuratedOvarianData.
# Data set from Tothill et al., 2008. Gene expression and clinical measurements with survival data.
# Models applied: AFT, Cox PH, Glmnet, RSF, GP, GPS1, GPS2, GPS3
#--------------------------------------------------------------------------------------------------------------------------------------------#

## RunAllModelsIncGPLaplace (2016_03_09) 																	##
## Applying Gaussian process to synthetic data generated with censored survival times in the training set. 	##
## Data provided. 10 repeats available, change nReps to run more. 											##
## All models applied to the same data. Models are:	GPNonSurvNoCens, GPSurvNoCorr, GPSurvCorrL, 			##
## 													GPSurvCorrV, GPSurvLaplace, AFT, Coxph, Glmnet, GBM, 	##
##													RF and RFSurvival. 										##
## Results saved to folder 'Runs' within working directory. 												##
## GPSurv with Laplace approximation is implemented by ApplyGPLaplace, and other GP methods by ApplyGP.		##
## Gradients are not available for optimisation within GPLaplace, but are available for other GP methods.	##
## To implement profiling, uncomment lines bracketing model calls and analysis at the end of this file.		##


##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
library(curatedOvarianData)
library(fields)
library(gbm)
library(glmnet)
library(impute)
library(ipred)
library(MASS)
library(Matrix)
library(mice)
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
source('ApplyRFSurvival.R')
source('CalculateMetrics.R')
source('CensorData.R')
source('CovFunc.R')
source('DataExtractionCuratedOvarian.R')
source('ImputeMissingData.R')
source('GenerateData.R')
source('GetDataSetFeatures.R')
source('LogPriorX.R')
source('MakeSyntheticData.R')
source('MeanFunc.R')
source('NormaliseExpressionData.R')
source('PlotKaplanMeier.R')
source('PreLearnHyperparam.R')
source('PrintOptions.R')
source('PrintOptions2.R')
source('RemoveCensored.R')
source('RenameGeneFeaturesAccordingToDataSetCuratedOvarian.R')
source('SetParametersRealTothill.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
nReps 			<- 30

dataSetChoice 	<- 'OCGS'
switch(dataSetChoice,
	'OCGS' 			={geneSubsetFlag 		<- 'TaqMan' 
					  clinicalSubsetFlag 	<- 'None'
					  dimension 			<- 97},
	'OCGS+Clinical'	={geneSubsetFlag 		<- 'TaqMan' 
					  clinicalSubsetFlag 	<- 'Four'
					  dimension 			<- 100},
	'SRGS'			={geneSubsetFlag 		<- 'SRGS1'
					  clinicalSubsetFlag 	<- 'None'
					  dimension 			<- 84},
	'SRGS+Clinical'	={geneSubsetFlag 		<- 'SRGS1'
					  clinicalSubsetFlag 	<- 'Four'
					  dimension 			<- 87})

allParameterStructures 	<- SetParametersRealTothill(geneSubsetFlag=geneSubsetFlag,clinicalSubsetFlag=clinicalSubsetFlag,dimension=dimension,nReps=nReps)
plotSaveOptions 		<- allParameterStructures$plotSaveOptions
dataOptionsStructure 	<- allParameterStructures$dataOptionsStructure
parameterStructure 		<- allParameterStructures$parameterStructure

censoringType 			<- dataOptionsStructure$censoringType

##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureGPNonSurvNoCens 	<- list()
outputStructureGPSurvNoCorr 	<- list()
outputStructureGPSurvCorrV		<- list()
outputStructureGPSurvCorrL 		<- list()
outputStructureAFT				<- list()
outputStructureCoxph 			<- list()
outputStructureGlmnet 			<- list()
outputStructureRFSurvival 		<- list()

c.index.GPNonSurvNoCens 		<- rep(NA,nReps)
c.index.GPSurvNoCorr 			<- rep(NA,nReps)
c.index.GPSurvCorrV				<- rep(NA,nReps)
c.index.GPSurvCorrL 			<- rep(NA,nReps)
c.index.AFT						<- rep(NA,nReps)
c.index.Coxph 					<- rep(NA,nReps)
c.index.Glmnet 					<- rep(NA,nReps)
c.index.RFSurvival 				<- rep(NA,nReps)

rmse.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.GPSurvNoCorr 				<- rep(NA,nReps)
rmse.GPSurvCorrV				<- rep(NA,nReps)
rmse.GPSurvCorrL 				<- rep(NA,nReps)
rmse.AFT						<- rep(NA,nReps)
rmse.Coxph 						<- rep(NA,nReps)
rmse.Glmnet 					<- rep(NA,nReps)
rmse.RFSurvival 				<- rep(NA,nReps)


##-------------------------------------------------------------------------------------##
##----------------------------------- Generate Data -----------------------------------##
##-------------------------------------------------------------------------------------##
trainingTestStructureForNReps 	<- GenerateData(dataOptionsStructure,outerFolder,nReps)


##-------------------------------------------------------------------------------------##
##----------------------------------- Run GP Model ------------------------------------##
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
##---------------------------------- Run Glmnet Model ---------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGlmnet[[i]] 			<- ApplyGlmnet(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Glmnet[i] 					<- ifelse(length(outputStructureGlmnet[[i]]$c.index)!=0,outputStructureGlmnet[[i]]$c.index,NA)
	rmse.Glmnet[i] 						<- ifelse(length(outputStructureGlmnet[[i]]$rmse)!=0,outputStructureGlmnet[[i]]$rmse,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------- Run RFSurvival Model --------------------------------##
##-------------------------------------------------------------------------------------##
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
if(plotSaveOptions$printResults){
	cat('---------------------------------------',fill=TRUE)
	cat('Concordance Indices:',fill=TRUE)
	cat('GPNonSurvNoCens mean c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPSurvNoCorr mean c index =',paste0(round(c.index.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrV mean c index =',paste0(round(c.index.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrL mean c index =',paste0(round(c.index.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('AFT mean c index =',paste0(round(c.index.AFT,4),collapse=', '),fill=TRUE)
	cat('Coxph mean c index =',paste0(round(c.index.Coxph,4),collapse=', '),fill=TRUE)
	cat('Glmnet mean c index =',paste0(round(c.index.Glmnet,4),collapse=', '),fill=TRUE)
	cat('RFSurvival mean c index =',paste0(round(c.index.RFSurvival,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
	cat('RMSE:',fill=TRUE)
	cat('GP mean rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS1 mean rmse =',paste0(round(rmse.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPS2 mean rmse =',paste0(round(rmse.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('GPS3 mean rmse =',paste0(round(rmse.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
}

							##-------------------------------##
							##--- Plot Kaplan-Meier Plots ---##
							##-------------------------------##
if(plotSaveOptions$savePlots){
	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvNoCorr','/',parameterStructure$unid,'GPS1','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvNoCorr[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrL','/',parameterStructure$unid,'GPS2','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrL[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrV','/',parameterStructure$unid,'GPS3','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrV[[i]]$plot3)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','AFT','/',parameterStructure$unid,'AFT','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureAFT[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','Coxph','/',parameterStructure$unid,'Coxph','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureCoxph[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','Glmnet','/',parameterStructure$unid,'Glmnet','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGlmnet[[i]]$plotKM)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','RFSurvival','/',parameterStructure$unid,'RFSurvival','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureRFSurvival[[i]]$plotKM)
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot Measured/Predicted ---##
							##-------------------------------##
if(plotSaveOptions$savePlots){
	pdf(paste0('Runs/',parameterStructure$unid,'/','GPNonSurv','/',parameterStructure$unid,'GP','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPNonSurvNoCens[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvNoCorr','/',parameterStructure$unid,'GPS1','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvNoCorr[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrL','/',parameterStructure$unid,'GPS2','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrL[[i]]$plot2)
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrV','/',parameterStructure$unid,'GPS3','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			replayPlot(outputStructureGPSurvCorrV[[i]]$plot2)
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot GP Hyperparameters ---##
							##-------------------------------##
if(plotSaveOptions$savePlots){
	pdf(paste0('Runs/',parameterStructure$unid,'/','GPNonSurv','/',parameterStructure$unid,'GP','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(rbind(unlist(t(outputStructureGPNonSurvNoCens[[i]]$logHypChosen[1:3])),unlist(t(outputStructureGPNonSurvNoCens[[i]]$parameterStructure$logHypStart[1:3])))
					/unlist(outputStructureGPNonSurvNoCens[[i]]$dataOptionsStructure$logHypGenerate)[1:3],
					type='l',ylab ='Hyperparameters',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvNoCorr','/',parameterStructure$unid,'GPS1','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvNoCorr[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvNoCorr[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrL','/',parameterStructure$unid,'GPS2','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvCorrL[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvCorrL[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',parameterStructure$unid,'/','GPSurvCorrV','/',parameterStructure$unid,'GPS3','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			matplot(t(t(outputStructureGPSurvCorrV[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvCorrV[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters',xlab='Cycle number')
			legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
		}
	dev.off(pdf.output)

}


							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##

modelNames 		<- c('AFT','Coxph','Glmnet','RFSurvival','GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV') # Reordered version of modelsList
modelNamesPlot 	<- c('AFT','Cox PH','Glmnet','RF Survival','GP','GPS1','GPS2','GPS3')
toInvert 		<- c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
for(i in 1:length(modelNames)){
	c.index.mat[,i] 	<- get(paste0('c.index.',modelNames[i]))
	if(toInvert[i]) c.index.mat[,i] <- 1-c.index.mat[,i]
	c.index.mean[i] 	<- mean(c.index.mat[,i],na.rm=TRUE)
}

if(plotSaveOptions$savePlots){
	pdf(file=paste0(getwd(),'/','Runs','/',parameterStructure$unid,'/',parameterStructure$unid,'PlotCIndexAllModels.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		plot(sort(rep(1:length(modelNames),nReps)),c(c.index.mat),pch=20,col=add.alpha('cadetblue3',0.8),cex=0.7,xaxt='n',xlab='',ylab='Mean Concordance Index',ylim=c(0.4,1))
		points(1:length(modelNames),c.index.mean,pch=20,cex=1.2)
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		text(1:length(modelNames),rep(0.97,length(modelNames)),labels=round(c.index.mean,4),cex=0.7,pos=3)
		layout(1)
	dev.off(pdf.output)

	pdf(file=paste0(getwd(),'/','Runs','/',parameterStructure$unid,'/',parameterStructure$unid,'PlotCIndexAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(c.index.mat,ylim=c(0.4,1),names=modelNamesPlot,ylab='Concordance Index',las=2,col='cadetblue3')
		layout(1)
	dev.off(pdf.output)
}


##-------------------------------------------------------------------------------------##
##------------------------------------ Save Output ------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureAll <- list('outputStructureGPNonSurvNoCens'=outputStructureGPNonSurvNoCens,
							'outputStructureGPSurvNoCorr'=outputStructureGPSurvNoCorr,
							'outputStructureGPSurvCorrV'=outputStructureGPSurvCorrV,
							'outputStructureGPSurvCorrL'=outputStructureGPSurvCorrL,
							'outputStructureAFT'=outputStructureAFT,
							'outputStructureCoxph'=outputStructureCoxph,
							'outputStructureGlmnet'=outputStructureGlmnet,
							'outputStructureRFSurvival'=outputStructureRFSurvival)

save(list='outputStructureAll',file=paste0('Runs','/',parameterStructure$unid,'/',parameterStructure$unid,'_','outputStructureAll','_',geneSubsetFlag,'_',clinicalSubsetFlag,'_','Workspace.RData'))
