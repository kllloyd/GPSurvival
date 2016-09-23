MeanFunc <-function(meanHyp,x,meanFuncForm,n){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# Mean function, Linear or Zero mean
	#-------------------------------------------------------------------------------------------------------#

	meanHyp <- as.matrix(meanHyp,nrow=length(meanHyp))
	dimension <- dim(x)[2]
  	switch(meanFuncForm,
  			'Linear' = {A = x%*%meanHyp[1:dimension,,drop=FALSE]+meanHyp[(dimension+1),]*matrix(1,n,1)},
  			'Zero'   = {A = matrix(0,n,1)}
  	)
  	return(A)
}