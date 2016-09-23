ApplyCoxYuanEtAl <- function(trainingTestStructure){
    #-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
    # Applies Cox PH model as Yuan et al. (2014)                                     
    #-------------------------------------------------------------------------------------------------------#
	source('coxcvMolClinical.R')
	source('coxcvMol.R')
	source('coxcvClinical.R')

	cancer 			<- trainingTestStructure$cancer
	molPlatform 	<- trainingTestStructure$molPlatform
	clinicalFlag 	<- trainingTestStructure$clinicalFlag

	x.train 		<- trainingTestStructure$x.train
	x.test 			<- trainingTestStructure$x.test
	y.train 		<- trainingTestStructure$y.train
	y.test 			<- trainingTestStructure$y.test
	clinical.train 	<- trainingTestStructure$clinical.train
	clinical.test 	<- trainingTestStructure$clinical.test

    nSamples        <- trainingTestStructure$nSamples
    dimension       <- trainingTestStructure$dimension

	if(clinicalFlag&molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    	output.cox <- try(coxcvMolClinical(x.train, y.train, x.test, y.test, clinical.train, clinical.test))
    } else if(clinicalFlag&!(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein'))){
    	output.cox <- try(coxcvClinical(clinical.train, y.train, clinical.test, y.test))
    } else if(!clinicalFlag&molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    	output.cox <- try(coxcvMol(x.train, y.train, x.test, y.test))
    }

    if (class(output.cox)=="try-error"){
        result.cox   <- rep(NA, nSamples)
        c.index.cox  <- NA
    } else {
        result.cox   <- output.cox$cox.predict
        c.index.cox  <- output.cox$c.index.predict
    }

    toReturn <- list('predictions'=result.cox,'c.index'=c.index.cox)

    return(toReturn)
}