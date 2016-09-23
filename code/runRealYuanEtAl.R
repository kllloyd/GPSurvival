#--------------------------------------------------------------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
# Recreating analysis of TCGA data by Yuan et al., 2014
# Cox PH, RSF and feature selection written using code from Synapse, ID:syn1720423, as reference. Where possible, original code has been used.
#--------------------------------------------------------------------------------------------------------------------------------------------#
# runRealYuanEtAl applies models to TCGA data curated by Yuan et al., 2014 to predict survival.
# Cox PH and RSF applied as Yuan et al. 
# Also applied, GP using only uncensored samples, GPS1, GPS2 and GPS3.
# Synapse login required when prompted to access data.
# Cancers: KIRC, LUSC, GBM, OV
# Platforms: Clinical, SCNA, mRNA, miRNA, protein (methyl is also available)
#--------------------------------------------------------------------------------------------------------------------------------------------#

library(fields)
library(foreach)
library(glmnet)
library(impute)
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
library(synapseClient)
synapseLogin('k.lloyd.1@warwick.ac.uk', 'Synapse3182')
# synapseLogin()

source('add.alpha.R')
source('ApplyFeatureSelection.R')
source('AdjustTrainingSurvivalMeanVariance.R')
source('ApplyCoxYuanEtAl.R')
source('ApplyGP.R')
source('ApplyRSFYuanEtAl.R')
source('CalculateMetrics.R')
source('CensorData.R')
source('CovFunc.R')
source('GenerateDataYuanEtAl.R')
source('GetSynapseIDs.R')
source('LogPriorX.R')
source('MeanFunc.R')
source('NormaliseExpressionData.R')
source('PlotKaplanMeier.R')
source('PreLearnHyperparam.R')
source('PrintOptions2.R')
source('read_table.R')
source('RemoveCensored.R')
source('RunAllModelsAllCombinations.R')
source('SetParametersRealYuan.R')

cancersMolecular 	<- list('KIRC'=c('None','SCNA','mRNA','miRNA','protein','SCNA','mRNA','miRNA','protein'),
							'OV'=c('None','SCNA','mRNA','miRNA','protein','SCNA','mRNA','miRNA','protein'),
							'GBM'=c('None','SCNA','mRNA','miRNA','SCNA','mRNA','miRNA'),
							'LUSC'=c('None','SCNA','mRNA','miRNA','protein','SCNA','mRNA','miRNA','protein'))
cancersClinical 	<- 	list('KIRC'=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE),
							'OV'=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE),
							'GBM'=c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE),
							'LUSC'=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE))
nReps 				<- 100
modelNames 			<- c('Cox','RSF','GP','GPS1','GPS2','GPS3')
c.index 			<- matrix(NA,ncol=length(modelNames),nrow=nReps)
cancer 				<- 'KIRC'

# for(j in 1:length(cancersMolecular[[cancer]])){
for(j in 1:1){
	assign(paste0('outputStructure',cancer),RunAllModelsAllCombinations(cancer,cancersMolecular[[cancer]][j],cancersClinical[[cancer]][j],nReps))
	save(list=paste0('outputStructure',cancer),file=paste0('Runs','/',eval(parse(text=paste0('outputStructure',cancer,'$parameterStructure$unid'))),'/','outputStructure',cancer,'_',cancersMolecular[[cancer]][j],'_',cancersClinical[[cancer]][j],'_','Workspace.RData'))

	unid 				<- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[3],'[[',1,']]$parameterStructure$unid')))
	for(k in 1:nReps){
		c.index[k,1] <- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[1],'[[',k,']]$c.index')))
		c.index[k,2] <- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[2],'[[',k,']]$c.index')))
		c.index[k,3] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[3],'[[',k,']]$c.index')))
		c.index[k,4] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[4],'[[',k,']]$c.index')))
		c.index[k,5] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[5],'[[',k,']]$c.index')))
		c.index[k,6] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[6],'[[',k,']]$c.index')))
	}

	pdf(file=paste0('Runs',"/",unid,'/','PlotCIndex',cancer,'_',cancersMolecular[[cancer]][j],'_',cancersClinical[[cancer]][j],'.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		boxplot(c.index,ylim=c(0.3,1),names=modelNames,ylab='Concordance Index',col='cadetblue3',las=2)
	dev.off(pdf.output)

	assign(paste0('outputStructure',cancer),list())
}
