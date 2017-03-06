# Yuan Yuan 2012-12-13
##############################################################################
# Runs the cox proportional hazards model using nfolds cross validation and
# returns the predicted 
#
# nfolds = number of folds in the cross validation
# dir = directory where the scripts needed will be found (lasso.R)
#
# ** requires glmnet packages
##############################################################################
library(foreach)
# library(doMC)
# registerDoMC(6)

coxcvMol <- function(x.train,y.train, x.test, y.test, useLASSO=TRUE) {
	source("lassoSurv.R")
        source("cox_screen.R")
        # browser()
        
        # use univariate cox model for pre-selection, only keep the significant ones
        cols.include <- cox.screen(y.train, x.train, top=sum(y.train[,2]))
        if (length(cols.include)==0)
          {
            stop("No feature passed the univariate cox screen: exit.")
          }
        x.train <- x.train[,cols.include]
        x.test <- x.test[,cols.include]
        print(paste("After univariate cox screen, features remain:", length(cols.include)))

        if(useLASSO & length(cols.include) > 5) # further shrink by LASSO, if only a few features, no need to use LASSO, note 5 is a quite arbitrary setting
          {
            # do LASSO without cross validation to get the features to include in the model
            # change x.train and x.test to only include those features
            cols.include <- c()
            iter <- 0
            while(length(cols.include)<1)
              {
                set.seed(iter+1)
                cols.include <- try(lassoSurv(x=x.train,y=y.train,above=0))
                if (class(cols.include)=="try-error")
                  {
                    print("Errors occur while calculating by LASSO, recalculating...")
                    cols.include <- c()
                  }
                iter <- iter+1
                if(iter> 1)
                  {
                    print(paste(length(cols.include)," features selected. Recalculated by LASSO:", iter))
                  }
                if(iter>100) # maximum number of iterations allowed
                  {
                    stop("No significant features can be selected by LASSO: exit.")
                  }
              }    
            print(paste("After LASSO, features remain:", length(cols.include)))
            x.train <- x.train[,cols.include]
            x.test <- x.test[,cols.include]
          }
        # print(paste("seed =", seed, ": final features:"))
        print(paste("final features:"))
        # print("------------------------------------------------")
        print(colnames(x.train))
        print("------------------------------------------------")

        # cox model for prediction
        if (length(cols.include)==1)
          {
            cox <- coxph(y.train~x.train)
          }else
        {
          cox <- coxph(y.train~., data= data.frame(x.train))
        }
        
        cox$coefficients[is.na(cox$coefficients)]=0 # convert NAs to zero if any       
        cox.predict <- as.matrix(x.test)%*%cox$coefficients

        c.index.predict <- concordance.index(cox.predict, y.test[,1], y.test[,2])$c.index
        print(c.index.predict)

        toReturn <- list('cox.predict'=cox.predict,'c.index.predict'=c.index.predict)
        
	# return(cox.predict)
      return(toReturn)
}
