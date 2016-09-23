ApplyRSFYuanEtAl <- function(trainingTestStructure){
    #-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
    # Applies random survival forest as Yuan et al. (2015)                                     
    #-------------------------------------------------------------------------------------------------------#
	source('rsfcvMolClinical.R')
	source('rsfcvMol.R')
	source('rsfcvClinical.R')

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
    	output.rsf <- try(rsfcvMolClinical(x.train, y.train, x.test, y.test, clinical.train, clinical.test))
    } else if(clinicalFlag&!(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein'))){
    	output.rsf <- try(rsfcvClinical(clinical.train, y.train, clinical.test, y.test))
    } else if(!clinicalFlag&molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
    	output.rsf <- try(rsfcvMol(x.train, y.train, x.test, y.test))
    }

    if (class(output.rsf)=="try-error"){
        result.rsf   <- rep(NA, nSamples)
        c.index.rsf  <- NA
    } else {
        result.rsf   <- output.rsf$rsf.predict
        c.index.rsf  <- output.rsf$c.index.predict
    }

    toReturn <- list('predictions'=result.rsf,'c.index'=c.index.rsf)

    return(toReturn)
}