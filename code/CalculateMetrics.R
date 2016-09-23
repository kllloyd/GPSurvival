CalculateMetrics <- function(predicted,targets,events){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Function calculates concordance index and rmse of predictions, given target values and events
    #-------------------------------------------------------------------------------------------------------#
	
	nTest 		<- length(predicted)
	tiesFlag 	<- any(duplicated(predicted))
	if(tiesFlag){
		tied 			<- duplicated(predicted)|rev(duplicated(rev(predicted)))
		predicted[tied] <- predicted[tied]+rnorm(sum(tied),mean=0,sd=10^-6)
	}

	cIndexStructure <- concordance.index(x=predicted,surv.time=targets,surv.event=events)
	c.index 		<- cIndexStructure$c.index
	rmse 			<- sqrt(sum((targets-predicted)^2)/nTest)

	toReturn 		<- list('c.index'=c.index,'rmse'=rmse)

	return(toReturn)
}