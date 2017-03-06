CalculateMetrics <- function(testPredictions,testTargets,testEvents,trainingPredictions,trainingTargets,trainingEvents){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Function calculates concordance index and rmse of predictions, given target values and testEvents
    #-------------------------------------------------------------------------------------------------------#
	
	library(survcomp)
    library(survAUC)

	nTest 			<- length(testPredictions)
	tiesFlag 		<- any(duplicated(testPredictions))
	if(tiesFlag){
		tied 					<- duplicated(testPredictions)|rev(duplicated(rev(testPredictions)))
		testPredictions[tied] 	<- testPredictions[tied]+rnorm(sum(tied),mean=0,sd=10^-6)
	}

	cIndexStructure <- concordance.index(x=testPredictions,surv.time=testTargets,surv.event=testEvents)
	c.index 		<- cIndexStructure$c.index
	
	rmse 			<- sqrt(sum((testTargets-testPredictions)^2)/nTest)

	brierStructure 	<- predErr(Surv(trainingTargets,trainingEvents),Surv(testTargets,testEvents),trainingPredictions,testPredictions,seq(min(trainingTargets),max(trainingTargets),length.out=200),type='brier',int.type='weighted')
	brier 			<- brierStructure$ierror

	toReturn 		<- list('c.index'=c.index,'rmse'=rmse,'brier'=brier)

	return(toReturn)
}