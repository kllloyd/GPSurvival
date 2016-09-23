CovFunc <-function(x1,x2,maternParam,varFuncSqHyp,lengthHyp,covFuncForm){
  #--------------------------------------------------------------------------#
  # K Lloyd 2016_09_16
  #--------------------------------------------------------------------------#
  # Inputs x1 and x2 are always [dimension,nSamples] or [dimension,1]
  #--------------------------------------------------------------------------#

if(identical(x1,x2)){
  distance <- rdist(t(x1)/as.numeric(lengthHyp))
} else {
  distance <- rdist(t(x1)/as.numeric(lengthHyp), t(x2)/as.numeric(lengthHyp))
}

  maternParamFlag = as.character(maternParam)
  switch(covFuncForm,
         'Matern' = {switch(maternParamFlag,
                           '1' = {covariance = as.numeric(varFuncSqHyp)*exp(-distance)},
                           '3' = {covariance = as.numeric(varFuncSqHyp)*(1+(sqrt(3)*distance))*exp(-(sqrt(3)*distance))},
                           '5' = {covariance = as.numeric(varFuncSqHyp)*(1+(sqrt(5)*distance)+(5*distance^2)/3)*exp(-(sqrt(5)*distance))}
                           )},
         'SqExp'  = {covariance = as.numeric(varFuncSqHyp)*exp(-((distance)^2)/2)},
         'ARD'    = {covariance = as.numeric(varFuncSqHyp)*exp(-((distance)^2)/2)})

  return(covariance)
}