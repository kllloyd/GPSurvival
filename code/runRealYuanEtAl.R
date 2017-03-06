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


##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##

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
synapseLogin()


##-------------------------------------------------------------------------------------##
##--------------------------------- Source Functions ----------------------------------##
##-------------------------------------------------------------------------------------##

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


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##

cancer 				<- 'KIRC'
nReps 				<- 100

cancersMolecular 	<- list('KIRC'	=c('None','SCNA','methyl','mRNA','miRNA','protein','SCNA','methyl','mRNA','miRNA','protein'),
							'OV'	=c('None','SCNA','methyl','mRNA','miRNA','protein','SCNA','methyl','mRNA','miRNA','protein'),
							'GBM'	=c('None','SCNA','methyl','mRNA','miRNA','SCNA','methyl','mRNA','miRNA'),
							'LUSC'	=c('None','SCNA','mRNA','miRNA','protein','SCNA','mRNA','miRNA','protein'))
cancersClinical 	<- list('KIRC'	=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
							'OV'	=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
							'GBM'	=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE),
							'LUSC'	=c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE))
k 					<- switch(cancer,
							'KIRC' 	= 1,
							'OV' 	= 2,
							'GBM' 	= 3,
							'LUSC' 	= 4)
modelNames 			<- c('Cox','RSF','GP','GPS1','GPS2','GPS3')

##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##

c.index 			<- matrix(NA,ncol=length(modelNames),nrow=nReps)

for(k in 1:length(cancersMolecular[[cancer]])){
	## Run models & save results ##
	assign(paste0('outputStructure',cancer),RunAllModelsAllCombinations(cancer,cancersMolecular[[cancer]][k],cancersClinical[[cancer]][k],nReps,k))
	save(list=paste0('outputStructure',cancer),file=paste0('Runs','/',eval(parse(text=paste0('outputStructure',cancer,'$parameterStructure$unid'))),'/','outputStructure',cancer,'_',cancersMolecular[[cancer]][k],'_',cancersClinical[[cancer]][k],'_','Workspace.RData'))

	# Extract c index #
	unid 				<- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[3],'[[',1,']]$parameterStructure$unid')))
	for(i in 1:nReps){
		c.index[i,1] <- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[1],'[[',i,']]$c.index')))
		c.index[i,2] <- eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[2],'[[',i,']]$c.index')))
		c.index[i,3] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[3],'[[',i,']]$c.index')))
		c.index[i,4] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[4],'[[',i,']]$c.index')))
		c.index[i,5] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[5],'[[',i,']]$c.index')))
		c.index[i,6] <- 1-eval(parse(text=paste0('outputStructure',cancer,'$outputStructure',modelNames[6],'[[',i,']]$c.index')))
	}

	# Save plots #
	pdf(file=paste0('Runs',"/",unid,'/','PlotCIndex',cancer,'_',cancersMolecular[[cancer]][k],'_',cancersClinical[[cancer]][k],'.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		boxplot(c.index,ylim=c(0.3,1),names=modelNames,ylab='Concordance Index',col='cadetblue3',las=2)
	dev.off(pdf.output)

	assign(paste0('outputStructure',cancer),list())
}
