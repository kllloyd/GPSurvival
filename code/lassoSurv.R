# Yuan Yuan 2012-12-13

### get columns to include based on LASSO without cross validation.
### ** requires the glmnet package

lassoSurv <- function(x,y,above=0) {
	require(glmnet)
	# get rid of NAs
	# keep = indices in x that are not NA
         
	for (i in 1:ncol(x)) {
		keep <- c(1:nrow(x))[!is.na(x[,i])]
		x <- x[keep,]
		#convert the columns to numeric
		x[,i] <- as.numeric(x[,i])
		y <- y[keep]
	}
	
	x <- as.matrix(x)
	fit.cv <- cv.glmnet(x=x,y=y,family="cox",alpha=1,standardize=FALSE,nfolds=5) #data have been already standardized 
        
	lambda <- fit.cv$lambda.min
	final <- glmnet(x=x,y=y,family="cox",alpha=1,lambda=lambda, standardize=FALSE) # glmnet is standardized by default, so turn it off
	coef.fit <- coef(final,s=lambda)[2:(ncol(x)+1)] 	# exclude intercept
	cols.include <- which(abs(coef.fit) > above)
	return(cols.include)
}
