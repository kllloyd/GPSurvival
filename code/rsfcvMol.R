# Yuan Yuan 2012-12-13
##############################################################################
# Runs the cox proportional hazards model using nfolds cross validation and
# returns the predicted risk score
#
# nfolds = number of folds in the cross validation
# dir = directory where the scripts needed will be found (lasso.R)
#
# ** requires glmnet packages
##############################################################################
library(foreach)
# library(doMC)
library(randomForestSRC)
# registerDoMC(6)

rsfcvMol <- function(x.train,y.train, x.test, y.test, ntree=1000) {
	source("lassoSurv.R")
        source("cox_screen.R")
        pred.rsf <- list()

        cols.include <- cox.screen(y.train, x.train, top=sum(y.train[,2]))
        if (length(cols.include)==0)
          {
            stop("No feature passed the univariate cox screen: exit.")
          }
        x.train <- x.train[,cols.include]
        x.test <- x.test[,cols.include]
        print(paste("After univariate cox screen, features remain:", length(cols.include)))
        
        
        #Random survival forest
        data <- data.frame(cbind(y.train[,1],y.train[,2], x.train))
        if(length(cols.include)==1)
          {
            x.test <- data.frame(x.test)
            colnames(x.test) <- "x.train" # so that the formula is compatible
          }
        colnames(data)[1]="time"
        colnames(data)[2]="status"
        
        
        # rf <- rsf(Surv(time, status)~., data=data, ntree=ntree, seed=-1)
        rf.all <- rfsrc(Surv(time,status)~.,data=data,ntree=ntree,seed=-1,importance=TRUE)
        feature.imp <- names(rf.all$importance)[which(rf.all$importance>0)] # the features with non-zero importance
        print(paste("After random forest, features remain:", length(feature.imp)))
        print(feature.imp)
        print("------------------------------------------------")

        # prediction by rsf
        # rsf.predict <- predict(rf, data.frame(x.test), seed=-1)$mortality
        # rsf.train <- predict(rf, data.frame(x.train), seed=-1)$mortality 
        rsf.both.train <- predict(rf.all, data.frame(x.train), seed=-1)$predicted
        rsf.both.predict <- predict(rf.all, data.frame(x.test), seed=-1)$predicted

        library(survcomp)
        c.index.train <- concordance.index(rsf.both.train, y.train[,1], y.train[,2])$c.index
        c.index.predict <- concordance.index(rsf.both.predict, y.test[,1], y.test[,2])$c.index
        # c.index.train <- estC(scoreY=rsf.both.train, timeX=y.train[,1], statusX=y.train[,2])
        # c.index.predict <- estC(scoreY=rsf.both.predict, timeX=y.test[,1], statusX=y.test[,2])
        print(c.index.predict)

    toReturn <- list('rsf.predict'=rsf.both.predict,'c.index.predict'=c.index.predict)
    return(toReturn)
	# return(rsf.predict)
}
