GenerateDataYuanEtAl <- function(dataOptionsStructure){
	#----------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#----------------------------------------------------------------------------------------------------------#
	# Extracts data from Synapse
	#----------------------------------------------------------------------------------------------------------#

	nReps 			<- dataOptionsStructure$nReps

	cancer 			<- dataOptionsStructure$cancer		# c('KIRC','OV','GBM','LUSC')
	molPlatform 	<- dataOptionsStructure$molPlatform   # c('SCNA','methyl','mRNA','miRNA','protein','None')
	clinicalFlag 	<- dataOptionsStructure$clinicalFlag	# c(TRUE,FALSE)

	exp.ID 			<- dataOptionsStructure$exp.ID
	clinical.ID 	<- dataOptionsStructure$clinical.ID
	surv.ID 		<- dataOptionsStructure$surv.ID
	train.ID 		<- dataOptionsStructure$train.ID
	test.ID 		<- dataOptionsStructure$test.ID

	surv.data 		<- myRead(surv.ID)

	if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
		exp.data <- myRead(exp.ID)
	} else {
		exp.data <- NULL
	}

	if(clinicalFlag){
		clinical.data  <- myRead.simple(clinical.ID)
		if('gender'%in%colnames(clinical.data)){
			clinical.data$gender 	<- ifelse(clinical.data$gender=="FEMALE",1, 0)
			clinical.data$gender	<- as.numeric(clinical.data$gender)
		}
		if('grade'%in%colnames(clinical.data)){
			clinical.data$grade 							<- as.character(clinical.data$grade)
			clinical.data$grade[clinical.data$grade=='GB'] 	<- 0.5
			clinical.data$grade[clinical.data$grade=='G1'] 	<- 1
			clinical.data$grade[clinical.data$grade=='G2'] 	<- 2
			clinical.data$grade[clinical.data$grade=='G3'] 	<- 3
			clinical.data$grade[clinical.data$grade=='G4'] 	<- 4
			clinical.data$grade[clinical.data$grade=='GX'] 	<- 2.5
			clinical.data$grade								<- as.numeric(clinical.data$grade)
		}
		if('stage'%in%colnames(clinical.data)&cancer=='KIRC'){
			clinical.data$stage 									<- as.character(clinical.data$stage)
			clinical.data$stage[clinical.data$stage=="Stage I"] 	<- 1
			clinical.data$stage[clinical.data$stage=="Stage II"]	<- 2
			clinical.data$stage[clinical.data$stage=="Stage III"]	<- 3
			clinical.data$stage[clinical.data$stage=="Stage IV"] 	<- 4
			clinical.data$stage 									<- as.numeric(clinical.data$stage)
		}
		if('stage'%in%colnames(clinical.data)&cancer=='OV'){
			clinical.data$stage 								<- as.character(clinical.data$stage)
			clinical.data$stage[clinical.data$stage=='IA'] 		<- 1.0
			clinical.data$stage[clinical.data$stage=='IB'] 		<- 1.3
			clinical.data$stage[clinical.data$stage=='IC'] 		<- 1.6
			clinical.data$stage[clinical.data$stage=='IIA'] 	<- 2.0
			clinical.data$stage[clinical.data$stage=='IIB'] 	<- 2.3
			clinical.data$stage[clinical.data$stage=='IIC'] 	<- 2.6
			clinical.data$stage[clinical.data$stage=='IIIA'] 	<- 3.0
			clinical.data$stage[clinical.data$stage=='IIIB'] 	<- 3.3
			clinical.data$stage[clinical.data$stage=='IIIC'] 	<- 3.6
			clinical.data$stage[clinical.data$stage=='IV'] 		<- 4
			clinical.data$stage 								<- as.numeric(clinical.data$stage)
		}
		if('stage'%in%colnames(clinical.data)&cancer=='LUSC'){
			clinical.data$stage 									<- as.character(clinical.data$stage)
			clinical.data$stage[clinical.data$stage=='Stage IA'] 	<- 1.0
			clinical.data$stage[clinical.data$stage=='Stage IB'] 	<- 1.3
			clinical.data$stage[clinical.data$stage=='Stage II'] 	<- 2.0
			clinical.data$stage[clinical.data$stage=='Stage IIA'] 	<- 2.0
			clinical.data$stage[clinical.data$stage=='Stage IIB'] 	<- 2.3
			clinical.data$stage[clinical.data$stage=='Stage IIIA'] 	<- 3.0
			clinical.data$stage[clinical.data$stage=='Stage IIIB'] 	<- 3.3
			clinical.data$stage 									<- as.numeric(clinical.data$stage)
		}
	} else {
		clinical.data <- NULL
	}

	train.all 	<- read.table(synapse.read(train.ID),header=F, stringsAsFactors=F)
	test.all  	<- read.table(synapse.read(test.ID),header=F, stringsAsFactors=F)

	mySurv      <- Surv(surv.data$OS_OS, surv.data$OS_vital_status,type='right')

	toReturn 	<- list()
	for(i in 1:nReps){
		train.samples   <- train.all[,i]
		test.samples    <- test.all[,i]
		if(clinicalFlag){
			sampleNames <- rownames(clinical.data)
		} else if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
			sampleNames <- rownames(exp.data)
		}
		train.row       <- match(train.samples, sampleNames)
		test.row        <- match(test.samples, sampleNames)

		y.train         <- mySurv[train.row]
		y.test          <- mySurv[test.row]
		if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')){
			x.train     <- exp.data[train.row, ]
			x.test      <- exp.data[test.row, ]
		} else {
			x.train 	<- NULL
			x.test 		<- NULL
		}
		if(clinicalFlag){
			clinical.train  <- clinical.data[train.row, ]
			clinical.test   <- clinical.data[test.row,]
		} else {
			clinical.train 	<- NULL
			clinical.test 	<- NULL
		}

		trainingTargets <- as.matrix(surv.data$OS_OS[train.row],ncol=1)
		testTargets 	<- as.matrix(surv.data$OS_OS[test.row],ncol=1)
		trainingEvents 	<- as.matrix(surv.data$OS_vital_status[train.row],ncol=1)
		testEvents 		<- as.matrix(surv.data$OS_vital_status[test.row],ncol=1)

		if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein') & clinicalFlag){
			trainingData 	<- cbind(clinical.train,x.train)
			testData 		<- cbind(clinical.test,x.test)
		} else if(!(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein')) & clinicalFlag){
			trainingData 	<- clinical.train
			testData 		<- clinical.test
		} else if(molPlatform%in%c('SCNA','methyl','mRNA','miRNA','protein') & !clinicalFlag){
			trainingData 	<- x.train
			testData 		<- x.test
		}

		dimension		<- dim(trainingData)[2]
		nTraining		<- dim(trainingData)[1]
		nTest 			<- dim(testData)[1]
		nSamples 		<- nTraining + nTest


		toReturn[[i]] 	<- list('nReps'=nReps,'cancer'=cancer,'molPlatform'=molPlatform,'clinicalFlag'=clinicalFlag,'x.train'=x.train,'y.train'=y.train,
								'x.test'=x.test,'y.test'=y.test,'clinical.train'=clinical.train,'clinical.test'=clinical.test,
								'trainingTargets'=trainingTargets,'trainingData'=trainingData,'events'=trainingEvents,'testTargets'=testTargets,
								'testData'=testData,'testEvents'=testEvents,'dimension'=dimension,'nSamples'=nSamples,'nTraining'=nTraining,'nTest'=nTest)
	}

	return(toReturn)
}