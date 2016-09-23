CensorData <- function(trainingTestStructure,dataOptionsStructure){
    #---------------------------------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #---------------------------------------------------------------------------------------------------------------------------------#
    # Function censors randomly selected target values from dataset produced by MakeSyntheticData.R
    # 3 types of censoring: NormalLoopSample, NormalLoopRerun and Normal. NormalLoopRerun is not non-informative.
    # Normal selects new censored values and if the new value is shorter than the known target it is replaced by the censored time (numCensored not always = nCensored)
    # NormalLoop carries out the same procedure but if censored value is not shorter a new value is chosen until it is (numCensored = ncensored)
    #---------------------------------------------------------------------------------------------------------------------------------#

    #----------------- Extract data and targets -----------------#
	x              = trainingTestStructure$data
	y              = trainingTestStructure$targets
    censoringSD    = dataOptionsStructure$censoringSD
    censoringMean  = dataOptionsStructure$censoringMean
    nCensored      = dataOptionsStructure$nCensored
    censoringType  = dataOptionsStructure$censoringType
    yTransformed   = y                                      # No transformation previously applied (log etc.)
	events         = numeric(length(y))+1
	numCensored    = 0

    #---------------------- Censor targets ----------------------#
    if(nCensored!=0){
	   switch(censoringType,
	   	'NormalLoopRerun'  = {indicesToCensor          = sample(1:length(y),nCensored)
	   					      events[indicesToCensor]  = 0
	   					      numCensored              = nCensored
	   					      for(i in 1:length(indicesToCensor)){
        				      	censoredValue    = abs(rnorm(1, mean = censoringMean, sd = censoringSD))
        				      	while(censoredValue>yTransformed[indicesToCensor[i]]){
        				      	  censoredValue  = abs(rnorm(1, mean = censoringMean, sd = censoringSD))
        				      	}
        				      	yTransformed[indicesToCensor[i]] = censoredValue
        				      }
                             },
        'NormalLoopSample' = {yToCensor         = yTransformed
                              indicesToCensor   = 1:length(yToCensor)
                              while(numCensored<nCensored){
                                toCensor        = sample(indicesToCensor,1)
                                censoredValue   = abs(rnorm(1, mean = censoringMean, sd = censoringSD))
                                if(yToCensor[toCensor]>censoredValue){
                                    yToCensor[toCensor] = censoredValue
                                    events[toCensor]    = 0
                                    indicesToCensor     = indicesToCensor[-which(indicesToCensor==toCensor)]
                                    numCensored         = numCensored+1
                                }
                              }
                              yTransformed      = yToCensor
                             },
        'Normal'	       = {indicesToCensor  = sample(1:length(y),nCensored)
	   		                  for(i in 1:length(indicesToCensor)){
        			             censoredValue = abs(rnorm(1, mean = censoringMean, sd = censoringSD))
        			             if(yTransformed[indicesToCensor[i]]>censoredValue){
        			                yTransformed[indicesToCensor[i]] = censoredValue
        			                events[indicesToCensor[i]]       = 0
        			                numCensored                      = numCensored + 1
        			             }
        			          }
                             },
        'None'             = {numCensored = 0}
        )
    }
    # Gaussian process works on -infinity to infinity scale so rescale back to log(time)
    y = yTransformed

    toReturn = list('data'=x,'targets'=y,'numCensored'=numCensored,'events'=events)

}