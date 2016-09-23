PlotKaplanMeier <- function(predictions,testTargets,testEvents,model){
	#-----------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-----------------------------------------------------------------------------#
    # Plot Kaplan Meier plot and apply log rank test 
    #-----------------------------------------------------------------------------#

	library('rms')
	library('survival')

	plotFlag <- 'both'

	if(model%in%c('GPSurvCorrL','GPSurvCorrV','GPSurvNoCorr','GPNonSurv','RF')) testTargets <- exp(testTargets)

	nTest 				<- length(predictions)
	survivalDF 			<- data.frame('predictions'=c(predictions),'survivalObj'=Surv(c(testTargets),c(testEvents)),stringsAsFactors=FALSE)
	survivalDFOrdered 	<- survivalDF[order(survivalDF$predictions),]

	survivalDFOrdered$group <- rep(1,nTest)
	survivalDFOrdered$group[round(nTest/2):nTest] <- 2
	survivalDFOrdered$group <- factor(survivalDFOrdered$group)

	if(plotFlag=='rms'){
		survivalCurve <- npsurv(survivalObj~group,data=survivalDFOrdered)
		testDifference <- survdiff(survivalObj~group,data=survivalDFOrdered,rho=0)
		survplot(survivalCurve,xlab="Time (days)",ylab="Proportion Survival",n.risk=TRUE,col.fill=c(add.alpha('cadetblue3',0.5),add.alpha('chartreuse3',0.5)))
		title(main=paste0('Chisq = ',round(testDifference$chisq,2),', p = ',format(pchisq(testDifference$chisq, length(testDifference$n)-1, lower.tail = FALSE),scientific=TRUE,digits=3)))
	} else if(plotFlag=='survival'){
		survivalCurve <- survfit(survivalObj~group,data=survivalDFOrdered)
		testDifference <- survdiff(survivalObj~group,data=survivalDFOrdered,rho=0)
		plot(survivalCurve,conf.int=TRUE,mark.time=TRUE,col=c('blue','green'),xlab="Time (days)",ylab="Proportion Survival",
			main=paste0('Chisq = ',round(testDifference$chisq,2),', p = ',format(pchisq(testDifference$chisq, length(testDifference$n)-1, lower.tail = FALSE),scientific=TRUE,digits=3)))
	} else if (plotFlag=='both'){
		survivalCurve <- npsurv(survivalObj~group,data=survivalDFOrdered)
		testDifference <- survdiff(survivalObj~group,data=survivalDFOrdered,rho=0)
		survplot(survivalCurve,xlab="Time (days)",ylab="Proportion Survival",n.risk=TRUE,col.fill=c(add.alpha('cadetblue3',0.5),add.alpha('chartreuse3',0.5)),
				col=c(add.alpha('cadetblue3',0.5),add.alpha('chartreuse3',0.5)),label.curves=FALSE)
		title(main=paste0('Chisq = ',round(testDifference$chisq,2),', p = ',format(pchisq(testDifference$chisq, length(testDifference$n)-1, lower.tail = FALSE),scientific=TRUE,digits=3)))
		lines(survivalCurve,mark.time=TRUE)
	}
}