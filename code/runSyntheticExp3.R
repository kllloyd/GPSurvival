#----------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# All models applied to the same data. Models are:	GP, GPS1, RF, RSF and Cox PH 
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Changing proportion of training set censored. 7 values.
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
debugSource('add.alpha.R')
debugSource('AdjustTrainingSurvivalMeanVariance.R')
debugSource('ApplyAFT.R')
debugSource('ApplyCoxph.R')
debugSource('ApplyGBM.R')
debugSource('ApplyGlmnet.R')
debugSource('ApplyGP.R')
debugSource('ApplyRF.R')
debugSource('ApplyRFSurvival.R')
debugSource('CalculateMetrics.R')
debugSource('CensorData.R')
debugSource('CovFunc.R')
debugSource('DataExtraction.R')
debugSource('GenerateData.R')
debugSource('ImputeMissingData.R')
debugSource('LogPriorX.R')
debugSource('MakeSyntheticData.R')
debugSource('MeanFunc.R')
debugSource('NormaliseExpressionData.R')
debugSource('PlotKaplanMeier.R')
debugSource('PreLearnHyperparam.R')
debugSource('PrintOptions.R')
debugSource('PrintOptions2.R')
debugSource('RemoveCensored.R')
debugSource('SetParametersExp3.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
nReps <- 30


##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureGPSurvNoCorr 	<- list()
outputStructureGPNonSurvNoCens 	<- list()
outputStructureRF 				<- list()
outputStructureRFSurvival 		<- list()
outputStructureCoxph			<- list()
c.index.GPSurvNoCorr 			<- rep(NA,nReps)
c.index.GPNonSurvNoCens 		<- rep(NA,nReps)
c.index.RF 						<- rep(NA,nReps)
c.index.RFSurvival 				<- rep(NA,nReps)
c.index.Coxph 					<- rep(NA,nReps)
rmse.GPSurvNoCorr 				<- rep(NA,nReps)
rmse.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.RF 						<- rep(NA,nReps)
rmse.RFSurvival 				<- rep(NA,nReps)
rmse.Coxph 						<- rep(NA,nReps)
trainingTestStructureForNReps 	<- list()

censoredProportionToRun 		<- c(0.01,0.10,0.30,0.50,0.70,0.90,0.98)
unids 							<- rep(NA,length(censoredProportionToRun))
count 							<- 1

for(k in 1:length(censoredProportionToRun)){
	## Generate parameters ##
	allParameterStructures 					<- SetParametersExp3(nTraining=500,nTest=100,hypGenerateNoise=0.20,censoredProportion=censoredProportionToRun[k],nReps=nReps)

	## Make data ##
	trainingTestStructureForNReps 			<- GenerateData(allParameterStructures$dataOptionsStructure,allParameterStructures$plotSaveOptions$outerFolder,nReps)
	for(i in 1:nReps){
	trainingTestStructureForNReps[[i]] 		<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
	cat('Data set for run',i,'normalised',fill=TRUE)
	}

	## Run GP ##
	allParameterStructures$dataOptionsStructure$censoringType 		<- 'None'
	allParameterStructures$parameterStructure$modelType 			<- 'non-survival'
	for(i in 1:nReps){
		trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
		trainingTestStructure 				<- RemoveCensored(trainingTestStructure,allParameterStructures$dataOptionsStructure)
		outputStructureGPNonSurvNoCens[[i]] <- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		c.index.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$c.index)!=0,outputStructureGPNonSurvNoCens[[i]]$c.index,NA)
		rmse.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$rmse)!=0,outputStructureGPNonSurvNoCens[[i]]$rmse,NA)
	}
	allParameterStructures$dataOptionsStructure$censoringType 		<- 'NormalLoopSample'
	allParameterStructures$parameterStructure$modelType 			<- 'survival'

	## Run GPS1 ##
	for(i in 1:nReps){
		trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
		outputStructureGPSurvNoCorr[[i]] 	<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		c.index.GPSurvNoCorr[i] 			<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$c.index)!=0,outputStructureGPSurvNoCorr[[i]]$c.index,NA)
		rmse.GPSurvNoCorr[i] 				<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$rmse)!=0,outputStructureGPSurvNoCorr[[i]]$rmse,NA)
	}

	## Run RF ##
	allParameterStructures$dataOptionsStructure$censoringType 		<- 'None'
	allParameterStructures$parameterStructure$modelType 			<- 'non-survival'
	for(i in 1:nReps){
		trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
		trainingTestStructure 				<- RemoveCensored(trainingTestStructure,allParameterStructures$dataOptionsStructure)
		outputStructureRF[[i]] 				<- ApplyRF(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		c.index.RF[i] 						<- ifelse(length(outputStructureRF[[i]]$c.index)!=0,outputStructureRF[[i]]$c.index,NA)
		rmse.RF[i] 							<- ifelse(length(outputStructureRF[[i]]$rmse)!=0,outputStructureRF[[i]]$rmse,NA)
	}
	allParameterStructures$dataOptionsStructure$censoringType 		<- 'NormalLoopSample'
	allParameterStructures$parameterStructure$modelType 			<- 'survival'

	## Run RSF ##
	for(i in 1:nReps){
		trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
		outputStructureRFSurvival[[i]] 		<- ApplyRFSurvival(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		c.index.RFSurvival[i] 				<- ifelse(length(outputStructureRFSurvival[[i]]$c.index)!=0,outputStructureRFSurvival[[i]]$c.index,NA)
		rmse.RFSurvival[i] 					<- ifelse(length(outputStructureRFSurvival[[i]]$rmse)!=0,outputStructureRFSurvival[[i]]$rmse,NA)
	}

	## Run Cox PH ##
	for(i in 1:nReps){
		trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
		outputStructureCoxph[[i]] 			<- ApplyCoxph(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		c.index.Coxph[i] 					<- ifelse(length(outputStructureCoxph[[i]]$c.index)!=0,outputStructureCoxph[[i]]$c.index,NA)
		rmse.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$rmse)!=0,outputStructureCoxph[[i]]$rmse,NA)
	}

	## Print c.index & rmse ##
	if(allParameterStructures$plotSaveOptions$printResults){
		cat('---------------------------------------',fill=TRUE)
		cat('GPS1 mean c index =',paste0(round(c.index.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
		cat('GPS1 mean rmse =',paste0(round(rmse.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
		cat('GP mean c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
		cat('GP mean rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
		cat('RF mean c index =',paste0(round(c.index.RF,4),collapse=', '),fill=TRUE)
		cat('RSF mean c index =',paste0(round(c.index.RFSurvival,4),collapse=', '),fill=TRUE)
		cat('Cox PH mean c index =',paste0(round(c.index.Coxph,4),collapse=', '),fill=TRUE)
		cat('---------------------------------------',fill=TRUE)
	}

	## Save plots ##
	if(allParameterStructures$plotSaveOptions$savePlots){
		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','GPSurvNoCorr','/',allParameterStructures$parameterStructure$unid,'GPS1','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureGPSurvNoCorr[[i]]$plot3)
			}
		dev.off(pdf.output)

		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','GPSurvNoCorr','/',allParameterStructures$parameterStructure$unid,'GPS1','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureGPSurvNoCorr[[i]]$plot2)
			}
		dev.off(pdf.output)

		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','GPSurvNoCorr','/',allParameterStructures$parameterStructure$unid,'GPS1','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				matplot(t(t(outputStructureGPSurvNoCorr[[i]]$logHypTable[,1:3])/unlist(outputStructureGPSurvNoCorr[[i]]$dataOptionsStructure$logHypGenerate)[1:3]),
						type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
				legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
			}
		dev.off(pdf.output)
	}

	if(allParameterStructures$plotSaveOptions$savePlots){
		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','GPNonSurv','/',allParameterStructures$parameterStructure$unid,'GP','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureGPNonSurvNoCens[[i]]$plot2)
			}
		dev.off(pdf.output)

		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','GPNonSurv','/',allParameterStructures$parameterStructure$unid,'GP','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				matplot(rbind(unlist(t(outputStructureGPNonSurvNoCens[[i]]$logHypChosen[1:3])),unlist(t(outputStructureGPNonSurvNoCens[[i]]$parameterStructure$logHypStart[1:3])))
						/unlist(outputStructureGPNonSurvNoCens[[i]]$dataOptionsStructure$logHypGenerate)[1:3],
						type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
				legend('topleft',legend=c(expression('log'*sigma[n]^2),expression('log'*sigma[f]^2),expression('log'*l)),col=1:3,pch=c(NA,NA,NA),lty=c(1,2,3))
			}
		dev.off(pdf.output)
	}

	if(allParameterStructures$plotSaveOptions$savePlots){
		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','RF','/',allParameterStructures$parameterStructure$unid,'RF','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureRF[[i]]$plotKM)
			}
		dev.off(pdf.output)
	}

	if(allParameterStructures$plotSaveOptions$savePlots){
		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','RFSurvival','/',allParameterStructures$parameterStructure$unid,'RFSurvival','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureRFSurvival[[i]]$plotKM)
			}
		dev.off(pdf.output)
	}

	if(allParameterStructures$plotSaveOptions$savePlots){
		pdf(paste0('Runs/',allParameterStructures$parameterStructure$unid,'/','Coxph','/',allParameterStructures$parameterStructure$unid,'Coxph','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
			for(i in c(1:nReps)){
				replayPlot(outputStructureCoxph[[i]]$plotKM)
			}
		dev.off(pdf.output)
	}

	## Rename and save outputStructures so not overwritten ##
	unids[count] 	<- allParameterStructures$parameterStructure$unid

	assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvNoCorr'),outputStructureGPSurvNoCorr)
	assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPNonSurvNoCens'),outputStructureGPNonSurvNoCens)
	assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureRF'),outputStructureRF)
	assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureRFSurvival'),outputStructureRFSurvival)
	assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureCoxph'),outputStructureCoxph)

	save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvNoCorr'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'GPSurvNoCorrWorkspace.RData'))
	save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPNonSurvNoCens'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'GPNonSurvNoCensWorkspace.RData'))
	save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureRF'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'RFWorkspace.RData'))
	save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureRFSurvival'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'RFSurvivalWorkspace.RData'))
	save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureCoxph'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'CoxphWorkspace.RData'))

	count 			<- count + 1
}


##-------------------------------------------------------------------------------------##
##----------------------------------- Plot Results ------------------------------------##
##-------------------------------------------------------------------------------------##
c.index.GPSurvNoCorr.all 		<- matrix(0,nrow=length(censoredProportionToRun),ncol=nReps)
c.index.GPNonSurvNoCens.all 	<- matrix(0,nrow=length(censoredProportionToRun),ncol=nReps)
c.index.RF.all 					<- matrix(0,nrow=length(censoredProportionToRun),ncol=nReps)
c.index.RFSurvival.all 			<- matrix(0,nrow=length(censoredProportionToRun),ncol=nReps)
c.index.Coxph.all 				<- matrix(0,nrow=length(censoredProportionToRun),ncol=nReps)

c.index.GPSurvNoCorr.mean 		<- numeric()
c.index.GPNonSurvNoCens.mean 	<- numeric()
c.index.RF.mean 				<- numeric()
c.index.RFSurvival.mean 		<- numeric()
c.index.Coxph.mean 				<- numeric()
c.index.GPSurvNoCorr.median 	<- numeric()
c.index.GPNonSurvNoCens.median 	<- numeric()
c.index.RF.median 				<- numeric()
c.index.RFSurvival.median 		<- numeric()
c.index.Coxph.median 			<- numeric()

modelNamesPlot 					<- c('Cox PH','RF','RF Survival','GP','GPS1')
toInvert 						<- c(FALSE,TRUE,FALSE,TRUE,TRUE)

for(i in 1:length(censoredProportionToRun)){
	for(j in 1:nReps){
		c.index.GPSurvNoCorr.all[i,j] 		<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvNoCorr[[',j,']]$c.index')))
		c.index.GPNonSurvNoCens.all[i,j] 	<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPNonSurvNoCens[[',j,']]$c.index')))
		c.index.RF.all[i,j] 				<- 1-eval(parse(text=paste0(unids[i],'outputStructureRF[[',j,']]$c.index')))
		c.index.RFSurvival.all[i,j] 		<- eval(parse(text=paste0(unids[i],'outputStructureRFSurvival[[',j,']]$c.index')))
		c.index.Coxph.all[i,j] 				<- eval(parse(text=paste0(unids[i],'outputStructureCoxph[[',j,']]$c.index')))
	}
	c.index.GPSurvNoCorr.mean[i] 			<- mean(c.index.GPSurvNoCorr.all[i,],na.rm=TRUE)
	c.index.GPNonSurvNoCens.mean[i] 		<- mean(c.index.GPNonSurvNoCens.all[i,],na.rm=TRUE)
	c.index.RF.mean[i] 						<- mean(c.index.RF.all[i,],na.rm=TRUE)
	c.index.RFSurvival.mean[i] 				<- mean(c.index.RFSurvival.all[i,],na.rm=TRUE)
	c.index.Coxph.mean[i] 					<- mean(c.index.Coxph.all[i,],na.rm=TRUE)
	c.index.GPSurvNoCorr.median[i] 			<- median(c.index.GPSurvNoCorr.all[i,],na.rm=TRUE)
	c.index.GPNonSurvNoCens.median[i] 		<- median(c.index.GPNonSurvNoCens.all[i,],na.rm=TRUE)
	c.index.RF.median[i] 					<- median(c.index.RF.all[i,],na.rm=TRUE)
	c.index.RFSurvival.median[i] 			<- median(c.index.RFSurvival.all[i,],na.rm=TRUE)
	c.index.Coxph.median[i] 				<- median(c.index.Coxph.all[i,],na.rm=TRUE)
}

forBoxplot2 							<- matrix(NA,nrow=30,ncol=7*5)
forBoxplot2[,c(1,6,11,16,21,26,31)] 	<- t(c.index.GPNonSurvNoCens.all)
forBoxplot2[,c(2,7,12,17,22,27,32)] 	<- t(c.index.GPSurvNoCorr.all)
forBoxplot2[,c(3,8,13,18,23,28,33)] 	<- t(c.index.RF.all)
forBoxplot2[,c(4,9,14,19,24,29,34)] 	<- t(c.index.RFSurvival.all)
forBoxplot2[,c(5,10,15,20,25,30,35)] 	<- t(c.index.Coxph.all)

positions <- c(1.5-3,1.5-1.5,1.5,1.5+1.5,1.5+3,15-3,15-1.5,15,15+1.5,15+3,40-3,40-1.5,40,40+1.5,40+3,75-3,75-1.5,75,75+1.5,75+3,105-3,105-1.5,105,105+1.5,105+3,135-3,135-1.5,135,135+1.5,135+3,147-3,147-1.5,147,147+1.5,147+3)

cols <- matrix(c('deepskyblue3','blue3','chartreuse3','darkgreen','white','deepskyblue3','blue3','chartreuse3','darkgreen','black'),nrow=5)
ltys <- c(1:length(censoredProportionToRun))

pdf(paste0('Runs/','PlotAllModelsCIndexBoxplotStage3ColChange.pdf'),width=12,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	opar <- par('mar')
	layout(rbind(1,2), heights=c(8,1))
	boxplot(forBoxplot2,at=positions,xaxt = "n",yaxt='n',
			col='white',medcol='white',whiskcol='white',staplecol='white',boxcol='white',outcol='white',las=1)
	axis(side=1,at=positions[c(3,8,13,18,23,28,33)],labels=censoredProportionToRun)
	lines(positions[c(1,6,11,16,21,26,31)],c.index.GPNonSurvNoCens.median,col=cols[1,2],lty=1)
	lines(positions[c(2,7,12,17,22,27,32)],c.index.GPSurvNoCorr.median,col=cols[2,2],lty=1,lwd=1.7)
	lines(positions[c(3,8,13,18,23,28,33)],c.index.RF.median,col=cols[3,2],lty=1)
	lines(positions[c(4,9,14,19,24,29,34)],c.index.RFSurvival.median,col=cols[4,2],lty=1,lwd=1.7)
	lines(positions[c(5,10,15,20,25,30,35)],c.index.Coxph.median,col=cols[5,2],lty=1,lwd=1.7)
	boxplot(forBoxplot2,at=positions,xaxt = "n",
			col=rep(cols[,1],7),outcol=rep(cols[,2],7),ylab='Concordance Index',xlab='Proportion of Training Set Censored',notch=FALSE,add=TRUE,outcex=0.85,las=1)
	par(mar=c(0,0,0,0))
	plot.new()
	legend('center',c('GP','GPS1','RF','RSF','Cox PH'),col=cols[,2],pch=c(15,15,15,15,15),lty=c(1,1,1,1,1),lwd=c(1,1.5,1,1.5,1.5),ncol=3,bty ="n")
	par(mar=opar)
	layout(1)
dev.off(pdf.output)

positions <- c(1.5-2.6,1.5-1.3,1.5,1.5+1.3,1.5+2.6,15-2.6,15-1.3,15,15+1.3,15+2.6,40-2.6,40-1.3,40,40+1.3,40+2.6,75-2.6,75-1.3,75,75+1.3,75+2.6,105-2.6,105-1.3,105,105+1.3,105+2.6,135-2.6,135-1.3,135,135+1.3,135+2.6,147-2.6,147-1.3,147,147+1.3,147+2.6)

pdf(paste0('Runs/','PlotAllModelsCIndexBoxplotStage3ColChangeNoWhiskers.pdf'),width=12,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	opar <- par('mar')
	layout(rbind(1,2), heights=c(8,1))
	boxplot(forBoxplot2,at=positions,xaxt = "n",yaxt='n',
			col='white',medcol='white',whiskcol='white',staplecol='white',boxcol='white',outcol='white',las=1)
	axis(side=1,at=positions[c(3,8,13,18,23,28,33)],labels=censoredProportionToRun)
	lines(positions[c(1,6,11,16,21,26,31)],c.index.GPNonSurvNoCens.median,col=cols[1,2],lty=1)
	lines(positions[c(2,7,12,17,22,27,32)],c.index.GPSurvNoCorr.median,col=cols[2,2],lty=1,lwd=1.7)
	lines(positions[c(3,8,13,18,23,28,33)],c.index.RF.median,col=cols[3,2],lty=1)
	lines(positions[c(4,9,14,19,24,29,34)],c.index.RFSurvival.median,col=cols[4,2],lty=1,lwd=1.7)
	lines(positions[c(5,10,15,20,25,30,35)],c.index.Coxph.median,col=cols[5,2],lty=1,lwd=1.7)
	boxplot(forBoxplot2,at=positions,xaxt = "n",
			col=rep(cols[,1],7),outcol=rep(cols[,2],7),ylab='Concordance Index',xlab='Proportion of Training Set Censored',notch=FALSE,add=TRUE,outcex=0.85,las=1,
			whisklty=0,staplelty=0,outline=FALSE)
	par(mar=c(0,0,0,0))
	plot.new()
	legend('center',c('GP','GPS1','RF','RSF','Cox PH'),col=cols[,2],pch=c(15,15,15,15,15),lty=c(1,1,1,1,1),lwd=c(1,1.5,1,1.5,1.5),ncol=3,bty ="n")
	par(mar=opar)
	layout(1)
dev.off(pdf.output)