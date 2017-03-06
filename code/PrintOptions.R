PrintOptions <- function(dataOptionsStructure,parameterStructure,modelsList){
	#-----------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-----------------------------------------------------------------------------#

	dataSource 			<- dataOptionsStructure$dataSource
	logHypGenerate 		<- dataOptionsStructure$logHypGenerate
	covFuncFormGen 		<- dataOptionsStructure$covFuncFormGen
	meanFuncFormGen 	<- dataOptionsStructure$meanFuncFormGen
	extraParamGen 		<- dataOptionsStructure$extraParamGen
	dimension 			<- dataOptionsStructure$dimension
	nSamples 			<- dataOptionsStructure$nSamples
	nTraining 			<- dataOptionsStructure$nTraining
	nTest 				<- dataOptionsStructure$nTest
	gridMinimum 		<- dataOptionsStructure$gridMinimum
	gridMaximum 		<- dataOptionsStructure$gridMaximum
	censoringType 		<- dataOptionsStructure$censoringType
	useARD 				<- dataOptionsStructure$useARD
	extraDimensions 	<- dataOptionsStructure$extraDimensions
	censoredProportion 	<- dataOptionsStructure$censoredProportion

	model 				<- parameterStructure$model
	meanFuncForm 		<- parameterStructure$meanFuncForm
	covFuncForm 		<- parameterStructure$covFuncForm
	extraParam 			<- parameterStructure$extraParam
	maxit 				<- parameterStructure$maxit
	optimType 			<- parameterStructure$optimType
	logHypStart 		<- parameterStructure$logHypStart
	unid 				<- parameterStructure$unid

	cat('-------------------------------------------------------',fill=TRUE)
	cat('Predicting y from x data by generating synthetic data.',fill=TRUE)
	cat('-------------------------------------------------------',fill=TRUE)
	cat('Data generation options:',fill=TRUE)
	cat('\tGeneration was using Gaussian processes',fill=TRUE)
	cat('\tData source was',dataSource,fill=TRUE)
	cat('\tSamples were generated with',dimension,'dimensions',fill=TRUE)
	cat(paste0('\t',nTraining),'training samples and',nTest,'test samples were used',fill=TRUE)
	if(covFuncFormGen=='ARD') cat(paste0('\t',extraDimensions),'dimensions of uninformative data were included',fill=TRUE)
	cat('\tGenerating log hyperparameters were',unlist(logHypGenerate),fill=TRUE)
	cat('\tSamples were produced on a grid between',gridMinimum,'and',gridMaximum,fill=TRUE)
	if(exists('censoringType')&exists('censoringLevel')) cat('\tCensoring was',censoringType,'with level sigma =',censoringLevel,fill=TRUE)
	if(exists('censoredProportion')) cat('\tCensoring proportion of training and test sets was',censoredProportion,fill=TRUE)
	if(covFuncFormGen=='Matern'){cat('\tCovariance function was',covFuncFormGen,'with Matern parameter set to',extraParamGen,fill=TRUE)
	} else {cat('\tCovariance function was',covFuncFormGen,fill=TRUE)}
	cat('\tMean function was',meanFuncFormGen,fill=TRUE)

	cat('Fitting options:',fill=TRUE)
	cat('\tModels applied are ',paste0(modelsList,collapse=', '),fill=TRUE)
	cat('Files associated with this run will be labelled with the UNID',unid,fill=TRUE)
	cat('-------------------------------------------------------',fill=TRUE)

}
