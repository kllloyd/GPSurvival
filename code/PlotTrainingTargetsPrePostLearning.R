PlotTrainingTargetsPrePostLearning <- function(outputStructure){
	#---------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #---------------------------------------------------------------------------------------#
    # Plots synthetic training set pre and post censoring and post learning
    #---------------------------------------------------------------------------------------#

	opar <- par('mar')
	par(mar=c(5.1,4.3,4.1,2.1))
	layout(rbind(1,2), heights=c(7,1))

	trainingData 					<- outputStructure$trainingTestStructure$trainingData
	trainingEvents 					<- outputStructure$trainingTestStructure$events
	trainingTargets 				<- outputStructure$trainingTestStructure$trainingTargets
	trainingTargetsPreCensoring 	<- outputStructure$trainingTestStructure$trainingTargetsPreCensoring
	trainingTargetsCensoredEachRun 	<- outputStructure$trainingTargetsCensoredEachRun
	trainingTargetsLearned 			<- outputStructure$trainingTestStructure$trainingTargetsLearned
	testData 						<- outputStructure$trainingTestStructure$testData
	xLine 							<- outputStructure$xLine
	logHypChosen 					<- outputStructure$logHypChosen
	funcMeanPred 					<- outputStructure$funcMeanPred
	varData 						<- outputStructure$varData

	meanLine 						<- list()
	varLine							<- list()
	fLine 							<- list()
	varLine2 						<- list()
	fLine2 							<- list()
	varLine3 						<- list()
	fLine3 							<- list()
	for(i in 1:length(outputStructure$modelPlotOutputStructure)){
		meanLine[[i]] 				<- outputStructure$modelPlotOutputStructure[[i]]$funcMeanPred
		varLine[[i]]				<- outputStructure$modelPlotOutputStructure[[i]]$varData
		varLine2[[i]]				<- outputStructure$modelPlotOutputStructure[[i]]$varFunction
		varLine3[[i]]				<- outputStructure$modelPlotOutputStructure[[i]]$varData + outputStructure$modelPlotOutputStructure[[i]]$varCorrection
		fLine[[i]] 					<- rbind(as.matrix(meanLine[[i]]+2*sqrt(varLine[[i]])),as.matrix(rev(meanLine[[i]]-2*sqrt(varLine[[i]]))))
		fLine2[[i]] 				<- rbind(as.matrix(meanLine[[i]]+2*sqrt(varLine2[[i]])),as.matrix(rev(meanLine[[i]]-2*sqrt(varLine2[[i]]))))
		if(!length(varLine3[[i]])==0) fLine3[[i]] <- rbind(as.matrix(meanLine[[i]]+2*sqrt(varLine3[[i]])),as.matrix(rev(meanLine[[i]]-2*sqrt(varLine3[[i]]))))
	}

	cols 							<- c('aliceblue','cadetblue4','chartreuse4','chartreuse3','steelblue3')
	pchs 							<- c(NA,15,16,16,18)
	cexs 							<- c(NA,1.3,1.4,1.4,1.4)

	plot(c(0,0),ylim=c(min(c(trainingTargets,trainingTargetsPreCensoring,trainingTargetsLearned,unlist(fLine))),max(c(trainingTargets,trainingTargetsPreCensoring,trainingTargetsLearned,unlist(fLine)))),xlim=c(min(trainingData,testData),max(trainingData,testData)),col='white',ylab='Targets',xlab='Data',main='',cex.lab=1.6,cex.axis=1.6,las=1)
	polygon(rbind(xLine,as.matrix(rev(xLine))),fLine[[(dim(trainingTargetsCensoredEachRun)[2]-1)]],col=add.alpha(cols[1],1))
	lines(xLine,meanLine[[(dim(trainingTargetsCensoredEachRun)[2]-1)]])
	points(trainingData[trainingEvents==0],trainingTargets[trainingEvents==0],pch=pchs[2],col=add.alpha(cols[2],0.7),cex=cexs[2])
	points(trainingData[trainingEvents==1],trainingTargets[trainingEvents==1],pch=pchs[3],col=add.alpha(cols[3],0.7),cex=cexs[3])
	points(trainingData[trainingEvents==0],trainingTargetsPreCensoring[trainingEvents==0],pch=pchs[4],col=add.alpha(cols[4],0.7),cex=cexs[4])
	points(trainingData[trainingEvents==0],trainingTargetsCensoredEachRun[,(dim(trainingTargetsCensoredEachRun)[2]-1)],pch=pchs[5],col=add.alpha(cols[5],0.7),cex=cexs[5])
	# points(testData,funcMeanPred,col='purple',pch=18)
	segments(x0=trainingData[trainingEvents==0],y0=trainingTargets[trainingEvents==0],x1=trainingData[trainingEvents==0],y1=apply(cbind(trainingTargetsCensoredEachRun[,(dim(trainingTargetsCensoredEachRun)[2]-1)],trainingTargetsPreCensoring[trainingEvents==0]),1,max),col=add.alpha('black',0.5))
	par(mar=c(0,0,0,0))
	plot.new()
	legend('center',c('Training, Censored, Pre-censoring','Training, Censored, Post-censoring','Training, Uncensored','Training, Censored, Learned'),col=c(add.alpha(cols[4],0.7),add.alpha(cols[c(2,3,5)],0.7)),pch=c(pchs[4],pchs[c(2,3,5)]),pt.cex=c(cexs[4],cexs[c(2,3,5)]),lty=c(0,0,0,0),ncol=2,bty ="n",cex=1.6)
	# legend('center',c('Training, Censored, Post-censoring','Training, Uncensored','Training, Censored, Learned'),col=add.alpha(cols[2:4],0.7),pch=c(20,20,20),pt.cex=c(1.6,1.6,1.6),lty=c(0,0,0),ncol=2,bty ="n",cex=1.6)
	par(mar=opar)

	layout(1)
}