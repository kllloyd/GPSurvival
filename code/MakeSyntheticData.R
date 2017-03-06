MakeSyntheticData <- function(dataOptionsStructure){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Function generates synthetic data of the chosen dimensions, using the mean and covariance functions and hyperparameters supplied
	# y targets produced are on the scale -infinity to infinity
	# If ARD is used, uninformative, randomly generated x data dimensions are added on after generating y targets
	#-------------------------------------------------------------------------------------------------------#

	library('Matrix')
	methodFlag <- 'likeRealData'

	#-------------------- Extract parameters --------------------#
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

	if(useARD){
		dimension 				<- dimension-extraDimensions
		logHypGenerate$length 	<- logHypGenerate$length[1:dimension]
	}

	#----------------------- Select x data -----------------------#
	switch(methodFlag,
	# METHOD TO MAKE DATA LIKE REAL-LIFE DATA
	# Trying to fulfil requirements to be like real data: each variable generated from normal distribution with mean=0, sd=1
	# Y values may need transforming to be on right scale
	'likeRealData' 	= {	nSamplesRounded <- nSamples
						xAll 			<- matrix(0,nSamplesRounded,nSamplesRounded)
						xAll 			<- sapply(1:dimension,function(x) rnorm(nSamplesRounded,mean=0,sd=1))
						K  				<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
				 		if(!all( eigen(K)$values >10^-13 )){
				   			cat('K is not positive-definite. Nearby p.d. matrix will be found.',fill=TRUE)
				   			nearPDStructure <- nearPD(K)
				   			if(!nearPDStructure$converged) {
				   			  cat('Warning: Procedure did not converge.',fill=TRUE)
				   			  browser()
				   			} else {
				   			  cat('Procedure converged.',fill=TRUE)
				   			  K <- as.matrix(nearPDStructure$mat)
				   			}
				   		}},
				   		
	# METHOD 1
	# Samples are generated on a grid of a certain size to cover the space & ensure two samples are not too close
	'method1' 		= {	nSamplesRounded 		<- ceiling(nSamples^(1/dimension))^dimension
						pointsInOneDimension 	<- seq(from=gridMinimum,to=gridMaximum,length.out=ceiling(nSamplesRounded^(1/dimension)))
						xAll 					<- as.matrix(expand.grid(rep(list(pointsInOneDimension), dimension)))
						K  						<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						if(!all( eigen(K)$values >0 )){
						  cat('K is not positive-definite. Nearby p.d. matrix will be found.',fill=TRUE)
						  nearPDStructure <- nearPD(K)
						  if(!nearPDStructure$converged) {
						    cat('Warning: Procedure did not converge.',fill=TRUE)
						  } else {
						    cat('Procedure converged.',fill=TRUE)
						    K <- as.matrix(nearPDStructure$mat)
						  }
						}},

	# METHOD 2
	# Samples are generated randomly until a set are chosen which produce a useable covariance matrix (if values are too close together K is not p.d.)
	'method2' 		= {	nSamplesRounded <- ceiling(nSamples^(1/dimension))^dimension
						K 				<- matrix(0,nSamplesRounded,nSamplesRounded)
						while(inherits(try(chol(K), silent=TRUE), "try-error")){
							xAll 	<- matrix(rnorm(nSamplesRounded*dimension),nrow=nSamplesRounded)
							K 		<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						}},

	# METHOD 3
	# Samples are generated as part of a loop, ensuring that each new sample does not stop the covariance matrix being p.d. THIS WILL TAKE A LONG TIME FOR LARGE nSamples!!
	'method3' 		= {	nSamplesRounded <- ceiling(nSamples^(1/dimension))^dimension
						nStart 			<- 26
						K = matrix(0,nStart,nStart)
						while(inherits(try(chol(K), silent=TRUE), "try-error")){
							xAll 	<- matrix(rnorm(nStart*dimension),nrow=nStart)
							K 		<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						}
						nIn = nStart
						while(nIn<nSamplesRounded){
							xAllTest 	<- rbind(xAll,as.matrix(rnorm(1*dimension),nrow=1))
							K 			<- CovFunc(xAllTest,xAllTest,xAllTest,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
							if(!inherits(try(chol(K), silent=TRUE), "try-error")){
								xAll 	<- rbind(xAll,as.matrix(rnorm(1*dimension),nrow=1))
								nIn 	<- nIn +1
								cat('Sample',nIn,'selected',fill=TRUE)
							}
						}},

	# METHOD 4
	# Samples are generated on a grid with noise to add scatter
	'method4' 		= {	nSamplesRounded 		<- ceiling(nSamples^(1/dimension))^dimension
						pointsInOneDimension 	<- seq(from=gridMinimum,to=gridMaximum,length.out=ceiling(nSamplesRounded^(1/dimension)))
						K <- matrix(0,nSamplesRounded,nSamplesRounded)
						while(inherits(try(chol(K), silent=TRUE), "try-error")){
							xAll 	<- as.matrix(expand.grid(rep(list(pointsInOneDimension), dimension)))+as.matrix(rnorm(nSamplesRounded*dimension),nrow=nSamplesRounded)
							K  		<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						}},

	# METHOD 5
	# Samples randomly generated, nearPD used to produce covariance matrix
	'method5' 		= {	nSamplesRounded <- ceiling(nSamples^(1/dimension))^dimension
						xAll 			<- matrix(rnorm(nSamplesRounded*dimension),nrow=nSamplesRounded)
						K  				<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						if(!all( eigen(K)$values >0 )){
						  cat('K is not positive-definite. Nearby p.d. matrix will be found.',fill=TRUE)
						  nearPDStructure <- nearPD(K)
						  if(!nearPDStructure$converged) {
						    cat('Warning: Procedure did not converge.',fill=TRUE)
						  } else {
						    cat('Procedure converged.',fill=TRUE)
						    K <- as.matrix(nearPDStructure$mat)
						  }
						}},

	# METHOD 6
	# Samples randomly generated from uniform distribution, nearPD used to produce covariance matrix
	'method6' 		= {	nSamplesRounded <- ceiling(nSamples^(1/dimension))^dimension
						xAll 			<- matrix(runif(nSamplesRounded*dimension,min=gridMinimum,max=gridMaximum),nrow=nSamplesRounded)
						K  				<- CovFunc(xAll,xAll,xAll,NA,extraParamGen,exp(logHypGenerate$func),exp(logHypGenerate$length),covFuncFormGen)
						if(!all( eigen(K)$values >10^-13 )){
						  cat('K is not positive-definite. Nearby p.d. matrix will be found.',fill=TRUE)
						  nearPDStructure <- nearPD(K)
						  if(!nearPDStructure$converged) {
						    cat('Warning: Procedure did not converge.',fill=TRUE)
						    browser()
						  } else {
						    cat('Procedure converged.',fill=TRUE)
						    K <- as.matrix(nearPDStructure$mat)
						  }
						}}
	)

	#--------------------- Generate y targets --------------------#
	mu 		<- MeanFunc(logHypGenerate$mean,xAll,meanFuncFormGen,nSamplesRounded)
	d1 		<- rnorm(nSamplesRounded)
	d2 		<- rnorm(nSamplesRounded)
	yAll 	<- mu + t(chol(K))%*%d1 + sqrt(exp(logHypGenerate$noise))*d2
	yAll 	<- exp(yAll)

	#--------------------- Randomise samples ---------------------#
	indices <- sample(1:(dim(yAll)[1]),nSamplesRounded)
	y  		<- yAll[indices,,drop=FALSE]
	x  		<- xAll[indices,,drop=FALSE]

	#-------------------- Add extra dimensions -------------------#

	if(useARD){
		extraXData 	<- matrix(runif(nSamples*extraDimensions,min=min(x),max=max(x)),nrow=nSamples)
		x 			<- cbind(x,extraXData)
	}

	toOutput <- list('data'=x,'targets'=y,'methodFlag'=methodFlag)
	return(toOutput)
}