AdjustTrainingSurvivalMeanVariance <- function(a,mu,sigma){
	#---------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#---------------------------------------------------------------------------------------------------------------------------------#
	# Function incorperates prior knowledge of training target censoring into the target value prediction.
	# Truncated PDF is approximated by normal PDF and used to estimate mean and variance of the distribution.
	# ONE SAMPLE #
	#---------------------------------------------------------------------------------------------------------------------------------#

	alpha 	= (a-mu)/sigma

	#---------------- Draw values for normal distribution using predictions of mean and variance from Gaussian process ---------------#
	survivalTime 			= seq(mu-50*sigma,mu+50*sigma,length=5000)
	oldPDF 					= dnorm(survivalTime,mean=mu, sd=sqrt(sigma))

	#----------------------------------------- Truncate distribution at known censored value -----------------------------------------#
	oldPDFtruncated 		= oldPDF[which(survivalTime>=a)]
	survivalTimeTruncated 	= survivalTime[(length(oldPDF)-length(oldPDFtruncated)+1):length(oldPDF)]
	
	#------ Calculate value of standard PDF and CDF evaluated at alpha = (truncation point - predicted mean)/predicted variance ------#
	phi_alpha 	= dnorm(alpha,mean=0,sd=1)
	Phi_alpha 	= pnorm(alpha,mean=0,sd=1)
	
	#----------------------------------- Calculate mean and variance of the truncated distribution -----------------------------------#
	exp_alpha 	= mu+sigma*(phi_alpha/(1-Phi_alpha))
	var_alpha	= sigma^2*(1-(phi_alpha/(1-Phi_alpha))*((phi_alpha/(1-Phi_alpha))-alpha))
  
	if((1-Phi_alpha) < 10^-12){
		exp_alpha 	= a
		var_alpha 	= 0
	} # Special case: if a is outside the non-zero (by numeric approximation) region of the PDF, CDF = 1 -> calculated mean and variance are NULL

  
	toReturn = list('mu'=exp_alpha,'sigma'=var_alpha)
	return(toReturn)

}
