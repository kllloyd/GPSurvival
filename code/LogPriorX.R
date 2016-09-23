LogPriorX <- function(hyp,x,data){
	#------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #------------------------------------------------------------------#
    # To apply gamma distributed priors to hyperparameters
    # Log
    #------------------------------------------------------------------#

	switch(hyp,
		'noise'	={k <- 2
				  t <- 1},
		'func'	={k <- 2
				  t <- 1},
		'length'={k <- 2
				  t <- diff(quantile(unlist(data),probs=c(0.05,0.95)))}
	)

	logPriorX <- -log(dgamma(1,k,t)) - k*log(t) + (k-1)*log(x) -x/t

	return(logPriorX)
}