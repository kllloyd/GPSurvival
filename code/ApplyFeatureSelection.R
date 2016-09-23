ApplyFeatureSelection <- function(trainingTestStructure,dataOptionsStructure,parameterStructure){
  #-------------------------------------------------------------------------------------------------------#
  # K Lloyd 2016_09_16
  #-------------------------------------------------------------------------------------------------------#
  # Applies feature selection via significance when run as features in Cox PH model (Yuan et al., 2014)
  # Also normalises/standardises training and test data                             
  # To be applied to data for GP models                                             
  #-------------------------------------------------------------------------------------------------------#

  logHypStart     <- parameterStructure$logHypStart

  y.train         <- trainingTestStructure$y.train
  y.test          <- trainingTestStructure$y.test
  x.train         <- trainingTestStructure$x.train
  x.test          <- trainingTestStructure$x.test
  clinical.train  <- trainingTestStructure$clinical.train
  clinical.test   <- trainingTestStructure$clinical.test

  clinicalFlag    <- dataOptionsStructure$clinicalFlag
  molPlatform     <- dataOptionsStructure$molPlatform

  if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    cols.include    <- cox.screen(y.train, x.train, top=sum(y.train[,2]))
    if (length(cols.include)==0){
      stop("No feature passed the univariate cox screen: exit.")
    }
    x.train <- x.train[,cols.include]
    x.test  <- x.test[,cols.include]
    cat("After univariate cox screen, features remain:", length(cols.include),fill=TRUE)
  }

  if(clinicalFlag&molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    trainingData  <- cbind(x.train, clinical.train)
    testData      <- cbind(x.test, clinical.test)
  } else if(!clinicalFlag&molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    trainingData  <- x.train
    testData      <- x.test
  } else if(clinicalFlag&!(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein'))){
    trainingData  <- clinical.train
    testData      <- clinical.test
  }

  dimension             <- dim(trainingData)[2] 

  toReturn              <- trainingTestStructure
  toReturn$trainingData <- trainingData
  toReturn$testData     <- testData
  toReturn$dimension    <- dimension

  toReturn              <- NormaliseExpressionData(toReturn,normaliseFlag=TRUE,winsoriseFlag=FALSE)

  toReturn$logHypStart  <- list('noise'=logHypStart$noise,'func'=logHypStart$func,'length'=logHypStart$length,'mean'=rep(logHypStart$mean[1],dimension+1))

  return(toReturn)
}