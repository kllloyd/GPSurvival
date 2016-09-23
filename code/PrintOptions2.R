PrintOptions2 <- function(dataOptionsStructure,parameterStructure){
    #-------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------#

	if(dataOptionsStructure$dataSource=='Generate'){
		cat('logHypGenerate =',unlist(dataOptionsStructure$logHypGenerate),fill=TRUE)
		cat('logHypStart =',unlist(parameterStructure$logHypStart),fill=TRUE)
		cat('censoringType =',dataOptionsStructure$censoringType,fill=TRUE)
		cat('censoringSD =',dataOptionsStructure$censoringSD,fill=TRUE)
		cat('censoringMean =',dataOptionsStructure$censoringMean,fill=TRUE)
		cat('censoredProportion =',dataOptionsStructure$censoredProportion,fill=TRUE)
	}
	if(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse') cat('dataSet =',dataOptionsStructure$dataSet,fill=TRUE)
	if(dataOptionsStructure$dataSource=='TCGAYuanEtAl') cat('cancer =',dataOptionsStructure$cancer,fill=TRUE)
	cat('dimension =',dataOptionsStructure$dimension,fill=TRUE)
	cat('nTraining =',dataOptionsStructure$nTraining,fill=TRUE)
	cat('nTest =',dataOptionsStructure$nTest,fill=TRUE)
	cat('tolerance =',parameterStructure$tolerance,fill=TRUE)
	if(dataOptionsStructure$dataSource=='CuratedOvarian'){
		cat('geneExpressionFlag =',dataOptionsStructure$geneExpressionFlag,fill=TRUE)
		if(dataOptionsStructure$geneExpressionFlag) cat('geneSubsetFlag =',dataOptionsStructure$geneSubsetFlag,fill=TRUE)
		cat('clinicalFeaturesFlag =',dataOptionsStructure$clinicalFeaturesFlag,fill=TRUE)
		if(dataOptionsStructure$clinicalFeaturesFlag[1]) cat('chosenClinicalFeatures =',dataOptionsStructure$chosenClinicalFeatures,fill=TRUE)
		if(dataOptionsStructure$clinicalFeaturesFlag[2]) cat('chosenExtraClinicalFeatures =',dataOptionsStructure$chosenExtraClinicalFeatures,fill=TRUE)
	}
	if(dataOptionsStructure$dataSource=='TCGASynapse'|dataOptionsStructure$dataSource=='TCGA2STAT'){
		cat('geneExpressionFlag =',dataOptionsStructure$geneExpressionFlag,fill=TRUE)
		if(dataOptionsStructure$geneExpressionFlag) cat('geneSubsetFlag =',dataOptionsStructure$geneSubsetFlag,fill=TRUE)
		cat('clinicalFeaturesFlag =',dataOptionsStructure$clinicalFeaturesFlag,fill=TRUE)
		if(dataOptionsStructure$clinicalFeaturesFlag) cat('chosenClinicalFeatures =',dataOptionsStructure$chosenClinicalFeatures,fill=TRUE)
	}
	if(dataOptionsStructure$dataSource=='TCGAYuanEtAl'){
		cat('molPlatform =',dataOptionsStructure$molPlatform,fill=TRUE)
		cat('clinicalFlag =',dataOptionsStructure$clinicalFlag,fill=TRUE)
	}
	if(dataOptionsStructure$dataSource=='CuratedOvarian'|dataOptionsStructure$dataSource=='TCGA2STAT'|dataOptionsStructure$dataSource=='TCGASynapse') cat('applyPCA =',dataOptionsStructure$applyPCA,fill=TRUE)
}