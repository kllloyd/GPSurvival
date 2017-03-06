CovFunc <-function(x1,x2,x,y,extraParam,varFuncSqHyp,lengthHyp,covFuncForm){
  #--------------------------------------------------------------------------#
  # K Lloyd 2016_09_16
  #--------------------------------------------------------------------------#
  # Inputs x1 and x2 are always [dimension,nSamples] or [dimension,1]
  # RF kernel implemented by Matthew Neal, accessed via package gaussianProcess
  # 
  # x = training data
  # y = training targets
  # x1 = training (KXX or KXX*) or test (KX*X*) data
  # x2 = training (KXX) or test (KXX* or KX*X*) data 
  #--------------------------------------------------------------------------#
  library(devtools)
  library(randomForest)
  ## UNCOMMENT THESE IF USING RF KERNEL ##
  # library(gaussianProcess)
  # library(cacheMan)
  # library(Rcpp)
  # library(RcppArmadillo)

  if(covFuncForm!='RF'){
    x1 <- t(x1)
    x2 <- t(x2)
    if(identical(x1,x2)){
      distance <- rdist(t(x1)/as.numeric(lengthHyp))
    } else {
      distance <- rdist(t(x1)/as.numeric(lengthHyp), t(x2)/as.numeric(lengthHyp))
    }
  }

  if(covFuncForm=='Matern') maternParamFlag <- as.character(extraParam)
  switch(covFuncForm,
         'Matern'       = {switch(maternParamFlag,
                                 '1' = {covariance <- as.numeric(varFuncSqHyp)*exp(-distance)},
                                 '3' = {covariance <- as.numeric(varFuncSqHyp)*(1+(sqrt(3)*distance))*exp(-(sqrt(3)*distance))},
                                 '5' = {covariance <- as.numeric(varFuncSqHyp)*(1+(sqrt(5)*distance)+(5*distance^2)/3)*exp(-(sqrt(5)*distance))}
                                 )},
         'SqExp'        = {covariance <- as.numeric(varFuncSqHyp)*exp(-((distance)^2)/2)},
         'ARD'          = {covariance <- as.numeric(varFuncSqHyp)*exp(-((distance)^2)/2)},
         'RF'           = {rf                     <- randomForest(x,y,num.trees=500)
                           rf.additional.params   <- create.rf.additional.params(rf)
                           rf.kernel              <- create.kernel.object("randomForest", additional_params=rf.additional.params)
                           if(identical(x1,x2)){
                            # x1=x2=trainingData or x1=x2=testData
                            covariance            <- get_covariance_matrix(kernel=rf.kernel, x=x1, sigma.n=0, hyper.params=numeric(0), cache=NULL)
                           } else{
                            # x1=trainingData & x2=testData
                            covariance            <- get_kstar_matrix(kernel=rf.kernel, data.to.predict=x2, training.data=x1, hyper.params= numeric(0))
                            covariance            <- t(covariance)
                           }}
        )

  return(covariance)
}
