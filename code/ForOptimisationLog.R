ForOptimisationLog <- function(logHyp,trainingData,trainingTargets,trainingEvents,meanFuncForm,dimension,maternParam,covFuncForm,nSamples,noiseCorr,logHypNoiseCorrection,imposePriors){
  #-------------------------------------------------------------------------------------------------------------------#
  # K Lloyd 2016_09_16
  #-------------------------------------------------------------------------------------------------------------------#
  # Function computes and returns log marginal likelihood for hyperparameter optimisation
  #-------------------------------------------------------------------------------------------------------------------#
  
  #--------------------------------------------- Extract hyperparameters ---------------------------------------------#
  logHyp                  <- matrix(logHyp,nrow=length(logHyp))
  if(noiseCorr=='noiseCorrLearned'){
    logHypNoiseCorrection <- logHyp[length(logHyp),,drop=FALSE]
    logHyp                <- logHyp[-length(logHyp),,drop=FALSE]
  }
  logHypNoise             <- logHyp[1,,drop=FALSE]
  logHypFunc              <- logHyp[2,,drop=FALSE]
  if(!(covFuncForm=='ARD')){logHypLength <- logHyp[3,,drop=FALSE]} else {logHypLength <- logHyp[3:(dimension+2),,drop=FALSE]}
  if(meanFuncForm=='Zero'){logHypMean <- rep(0,dimension+1)} else {logHypMean <- logHyp[(length(logHyp)-dimension):length(logHyp),,drop=FALSE]}

  #-------------- Compute mean vector, covariance matrix and log marginal likelihood using training data -------------#
  meanTraining      <- MeanFunc(logHypMean,trainingData,meanFuncForm,nSamples)
  K                 <- CovFunc(t(trainingData),t(trainingData),maternParam,exp(logHypFunc),exp(logHypLength),covFuncForm)+diag(as.numeric(exp(logHypNoise)),dim(trainingData)[1],dim(trainingData)[1])
  if(noiseCorr=='noiseCorrVec'|noiseCorr=='noiseCorrLearned'){
    varNoiseSqCorrection                    <- rep(0,dim(K)[1])
    varNoiseSqCorrection[!trainingEvents]   <- exp(logHypNoiseCorrection)
    K                                       <- K+diag(as.numeric(varNoiseSqCorrection),dim(K)[1],dim(K)[2])
  }
  if(class(try(chol(K),silent=TRUE))!='try-error'){
    L <- chol(K)
  } else if(class(try(chol(K),silent=TRUE))=='try-error' & class(try(nearPD(K),silent=TRUE))!='try-error'){
    nearMatrixStructure   <- nearPD(K)
    nearMatrix            <- as.matrix(nearMatrixStructure$mat)
    nearMatrixConvergence <- nearMatrixStructure$converged
    if(nearMatrixConvergence){
      L <- chol(nearMatrix)
      # cat('Nearby matrix found and used for cholesky decomposition of covariance matrix',fill=TRUE)
    } else {stop('Error: Covariance matrix is not positive definite. Generation of nearby p.d. matrix failed.')}
  } else {stop('Error: Covariance matrix is not positive definite. Generation of nearby p.d. matrix failed.')}
  # If covariance matrix + noise hyperparameter is not positive definite, a nearby p.d. matrix is found. If this is not found, an error message is returned.

  L         <- t(L)     # Making lower triangular, not upper triangular
  alpha     <- solve(t(L),(solve(L,(trainingTargets-meanTraining))))

  logPriorSigmaN  <- LogPriorX('noise',exp(logHypNoise),trainingData)
  logPriorSigmaF  <- LogPriorX('func',exp(logHypFunc),trainingData)
  logPriorL       <- LogPriorX('length',exp(logHypLength),trainingData)

  if(imposePriors){
    logMarLik <- (1/2)*t(trainingTargets-meanTraining)%*%alpha+sum(log(diag(L)))+(nSamples/2)*log(2*pi) - logPriorSigmaN - logPriorSigmaF - logPriorL
  } else {
    logMarLik <- (1/2)*t(trainingTargets-meanTraining)%*%alpha+sum(log(diag(L)))+(nSamples/2)*log(2*pi)
  }

  # cat('Parameters =',unlist(logHyp),fill=TRUE)
  return(logMarLik)
}
