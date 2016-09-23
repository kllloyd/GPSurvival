OptimisationGradientLog <- function(logHyp,trainingData,trainingTargets,trainingEvents,meanFuncForm,dimension,maternParam,covFuncForm,nSamples,noiseCorr,logHypNoiseCorrection) {
  #---------------------------------------------------------------------------------------------#
  # K Lloyd 2016_09_16
  #---------------------------------------------------------------------------------------------#
  # Function computes and returns gradients for hyperparameter optimisation
  #---------------------------------------------------------------------------------------------#

  #---------------------------------- Extract hyperparameters ----------------------------------#
  logHyp       = matrix(logHyp,nrow=length(logHyp))
  if(noiseCorr=='noiseCorrLearned'){
    logHypNoiseCorrection = as.numeric(logHyp[length(logHyp),,drop=FALSE])
    logHyp                = logHyp[-length(logHyp),,drop=FALSE]
  }
  logHypNoise   = as.numeric(logHyp[1,,drop=FALSE])
  logHypFunc    = as.numeric(logHyp[2,,drop=FALSE])
  if(!(covFuncForm=='ARD')){logHypLength = as.numeric(logHyp[3,,drop=FALSE])} else {logHypLength = as.numeric(logHyp[3:(dimension+2),,drop=FALSE])}
  if(meanFuncForm=='Zero'){logHypMean = rep(0,dimension+1)} else {logHypMean = as.numeric(logHyp[(length(logHyp)-dimension):length(logHyp),drop=FALSE])}

  #------- Compute required components for gradients using mean and covariance functions -------#
  meanTraining = MeanFunc(logHypMean,trainingData,meanFuncForm,nSamples)
  r			       = as.matrix(dist(trainingData))
  K            = CovFunc(t(trainingData),t(trainingData),maternParam,exp(logHypFunc),exp(logHypLength),covFuncForm)
  Id           = diag(1,dim(K)[1],dim(K)[2])
  if(noiseCorr=='noiseCorrLearned'|noiseCorr=='noiseCorrVec'){
    Ic                                    = diag(1-c(trainingEvents),dim(K)[1],dim(K)[2])
    varNoiseSqCorrection                  = rep(0,dim(K)[1])
    varNoiseSqCorrection[!trainingEvents] = exp(logHypNoiseCorrection)
    KyInverse                             = solve(K+Id*exp(logHypNoise)+diag(varNoiseSqCorrection,dim(K)[1],dim(K)[2]))
  } else {
    KyInverse                             = solve(K+Id*exp(logHypNoise))
  }

  #------------- Compute gradients, given the covariance and mean functions chosen -------------#
  switch(covFuncForm,
    'SqExp'   = {r              = as.matrix(dist(trainingData))
                 dvarNoiseSqHyp = (1/2*sum(diag(KyInverse%*%Id))-1/2*t(trainingTargets)%*%KyInverse%*%Id%*%KyInverse%*%trainingTargets)*exp(logHypNoise)
                 dvarFuncSqHyp  = (1/2*sum(diag(KyInverse%*%exp((-r^2)/(2*exp(logHypLength)^2))))-1/2*t(trainingTargets)%*%KyInverse%*%exp((-r^2)/(2*exp(logHypLength)^2))%*%KyInverse%*%trainingTargets)*exp(logHypFunc)
                 dl             = (1/2*sum(diag(KyInverse%*%(exp(logHypFunc)*(r^2/exp(logHypLength)^3)*exp((-r^2)/(2*exp(logHypLength)^2)))))-1/2*t(trainingTargets)%*%KyInverse%*%(exp(logHypFunc)*(r^2/exp(logHypLength)^3)*exp((-r^2)/(2*exp(logHypLength)^2)))%*%KyInverse%*%trainingTargets)*exp(logHypLength)
                 dmn1           = -(t(trainingTargets-meanTraining)%*%KyInverse%*%matrix(1,dim(trainingData)[1],1))
                 dmn            = matrix(0,dimension,1)
                 for(i in 1:dimension){dmn[i,] = -(t(trainingTargets-meanTraining)%*%KyInverse%*%trainingData[,i,drop=FALSE])}
                 if(noiseCorr=='noiseCorrLearned') dvarNoiseSqCorrection = (1/2*sum(diag(KyInverse%*%Ic))-1/2*t(trainingTargets)%*%KyInverse%*%Ic%*%KyInverse%*%trainingTargets)*exp(logHypNoiseCorr)},
    'ARD'     = {dvarNoiseSqHyp = (1/2*sum(diag(KyInverse%*%Id))-1/2*t(trainingTargets)%*%KyInverse%*%Id%*%KyInverse%*%trainingTargets)
                 dvarFuncSqHyp  = (1/2*sum(diag(KyInverse%*%exp((-r^2)/2)))-1/2*t(trainingTargets)%*%KyInverse%*%exp((-r^2)/2)%*%KyInverse%*%trainingTargets)
                 dl             = matrix(0,dimension,1)
                 dmn            = matrix(0,dimension,1)
                 dsigmadl       = list(length(logHypLength))
                 for(k in 1:length(logHypLength)){
                  # dsigmadl[[k]] = matrix(0,dim(trainingData)[1],dim(trainingData)[1])
                  # for(i in 1:dim(trainingData)[1]){
                  #   for(j in 1:dim(trainingData)[1]){
                  #     dsigmadl[[k]][i,j] = exp(logHypFunc)*exp((1/4*t(trainingData[i,]-trainingData[j,])%*%t(t(trainingData[i,]-trainingData[j,])))/((exp(logHypLength)[k])^3))
                  #     if(i==j) dsigmadl[[k]][i,j] = dsigmadl[[k]][i,j]+exp(logHypNoise)
                  #   }
                  # }
                  # if(i==j) dsigmadl[[k]][i,j] = dsigmadl[[k]][i,j]+exp(logHypNoise)
                  dsigmadl[[k]] = exp(logHypFunc)*exp((1/4*as.matrix(dist(trainingData)))/((exp(logHypLength)[k])^3))+diag(exp(logHypNoise),nrow=dim(trainingData)[1],ncol=dim(trainingData)[1])
                 }
                 for(i in 1:dimension){dmn[i,] = -(t(trainingTargets-meanTraining)%*%KyInverse%*%trainingData[,i,drop=FALSE])
                                       dl[i]   = (1/2*sum(diag(KyInverse%*%dsigmadl[[i]]))-1/2*t(trainingTargets)%*%KyInverse%*%dsigmadl[[i]]%*%KyInverse%*%trainingTargets)
                                      }
                 dmn1           = -(t(trainingTargets-meanTraining)%*%KyInverse%*%matrix(1,dim(trainingData)[1],1))}
  )

  #------------------------- Form vector of gradients for optimisation -------------------------#
  if(meanFuncForm=='Zero') gradients = rbind(dvarNoiseSqHyp,dvarFuncSqHyp,dl) else gradients = rbind(dvarNoiseSqHyp,dvarFuncSqHyp,dl,dmn,dmn1)
  if(noiseCorr=='noiseCorrLearned') gradients = rbind(gradients,dvarNoiseSqCorrection)

  cat('Gradients =',gradients,fill=TRUE)

  return(gradients)
}