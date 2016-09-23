DataExtractionCuratedOvarian <- function(dataSet='TCGA',chosenClinicalFeatures,chosenExtraClinicalFeatures,chosenExpressionFeatures,geneExpressionFlag,clinicalFeaturesFlag){
	#------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_09_16
	#------------------------------------------------------------------------------------------------#
	# Extracts molecular and clinical data from R package curatedOvarianData
	#------------------------------------------------------------------------------------------------#

	library(curatedOvarianData)
	normaliseGeneExpression <- FALSE
	normaliseTo 			<- 'HMBS'

	data(list=paste0(dataSet,'_eset'))
	clinicalData <- pData(get(paste0(dataSet,'_eset')))
	# write.table(clinicalData,file=paste0('clinicalData_',dataSet,'.csv'),quote=FALSE,sep='^')
	# clinicalData 	= read.table(paste0('clinicalData_',dataSet,'.csv'),header=TRUE,stringsAsFactors=FALSE,sep='^',na.strings='NA')

	if(geneExpressionFlag){
		expressionData <- exprs(get(paste0(dataSet,'_eset')))
		if(length(chosenExpressionFeatures)==1){
			if(chosenExpressionFeatures=='All'){chosenExpressionFeatures 	<- rownames(expressionData)
				
			} else {
				renamedGeneFeatures 		<- RenameGeneFeaturesAccordingToDataSetCuratedOvarian(expressionData,chosenExpressionFeatures)
				chosenExpressionFeatures 	<- renamedGeneFeatures$geneListPresentInDataSet
			}
		} else {
			renamedGeneFeatures 		<- RenameGeneFeaturesAccordingToDataSetCuratedOvarian(expressionData,chosenExpressionFeatures)
			chosenExpressionFeatures 	<- renamedGeneFeatures$geneListPresentInDataSet
		}
		expressionFeatures 	<- as.data.frame(expressionData,stringsAsFactors=FALSE)[chosenExpressionFeatures,,drop=FALSE]
		names(expressionFeatures)=gsub('///','_',names(expressionFeatures))

		if(normaliseGeneExpression){
			normaliseTo 	<- RenameGeneFeaturesAccordingToDataSetCuratedOvarian(expressionData,normaliseTo)$geneListPresentInDataSet
			if(!is.na(normaliseTo)){
				temp 				<- as.data.frame(sapply(1:dim(expressionFeatures)[2],function(x) expressionFeatures[,x]/expressionData[normaliseTo,x]),stringsAsFactors=FALSE)
				row.names(temp) 	<- row.names(expressionFeatures)
				names(temp) 		<- names(expressionFeatures)
				expressionFeatures 	<- temp
			} else { 
				cat('Gene chosen for normalisation is not present in data set.\n Normalisation not carried out.',fill=TRUE) # This currently isn't reached, need to catch error in function too
			}
		}
		expressionFeatures <- t(expressionFeatures)
	}

	if(clinicalFeaturesFlag[1]){
		clinicalFeatures <- clinicalData
		clinicalFeatures 																					<- clinicalFeatures[,chosenClinicalFeatures]
		names(clinicalFeatures)[which(names(clinicalFeatures)=='vital_status')] 							<- 'event'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='days_to_death')] 							<- 'days_to_event'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='age_at_initial_pathologic_diagnosis')] 		<- 'age'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='days_to_tumor_recurrence')] 				<- 'days_to_recurrence'

		if('tax'%in%chosenClinicalFeatures){
			names(clinicalFeatures)[which(names(clinicalFeatures)=='tax')] 									<- 'taxane'
			clinicalFeatures$taxane[which(clinicalFeatures$taxane=='y')] 									<- 1
			clinicalFeatures$taxane[which(clinicalFeatures$taxane=='n')] 									<- 0
			clinicalFeatures$taxane 																		<- as.numeric(clinicalFeatures$taxane)
		}
		if('pltx'%in%chosenClinicalFeatures){
			names(clinicalFeatures)[which(names(clinicalFeatures)=='pltx')] 								<- 'platinum'
			clinicalFeatures$platinum[which(clinicalFeatures$platinum=='y')] 								<- 1
			clinicalFeatures$platinum[which(clinicalFeatures$platinum=='n')] 								<- 0
			clinicalFeatures$platinum 																		<- as.numeric(clinicalFeatures$platinum)
		}
		if('debulking'%in%chosenClinicalFeatures){
			clinicalFeatures$debulking[which(clinicalFeatures$debulking=='optimal')] 						<- 1
			clinicalFeatures$debulking[which(clinicalFeatures$debulking=='suboptimal')] 					<- 0
			clinicalFeatures$debulking 																		<- as.numeric(clinicalFeatures$debulking)
		}

		if('summarygrade'%in%chosenClinicalFeatures){
			clinicalFeatures$summarygrade[which(clinicalFeatures$event=='high')]							<- 1
			clinicalFeatures$summarygrade[which(clinicalFeatures$event=='low')]								<- 0
			clinicalFeatures$summarygrade 																	<- as.numeric(clinicalFeatures$summarygrade)
		}

		if(all(c('tumorstage','substage')%in%chosenClinicalFeatures)){
			clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='a')] 								<- clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='a')]
			clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='b')] 								<- clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='b')] + 0.3
			clinicalFeatures$tumorstage[which(clinicalFeatures$substage==NA)] 								<- clinicalFeatures$tumorstage[which(clinicalFeatures$substage==NA)] + 0.3
			clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='c')] 								<- clinicalFeatures$tumorstage[which(clinicalFeatures$substage=='c')] + 0.6
		}

		if('recurrence_status'%in%chosenClinicalFeatures){
			clinicalFeatures$recurrence_status[which(clinicalFeatures$recurrence_status=='recurrence')] 	<- 1
			clinicalFeatures$recurrence_status[which(clinicalFeatures$recurrence_status=='norecurrence')] 	<- 0
			clinicalFeatures$recurrence_status 																<- as.numeric(clinicalFeatures$recurrence_status)
		}

		if(clinicalFeaturesFlag[2]){
			temp 				<- strsplit(clinicalData$uncurated_author_metadata,'///')
			tempMatrix 			<- matrix(0,length(temp),length(temp[[1]]))
			extraClinicalData 	<- matrix(0,length(temp),length(temp[[1]]))

			for(i in 1:length(temp)){
				tempMatrix[i,] <- temp[[i]]
			}

			for(i in 1:dim(tempMatrix)[1]){
				extraClinicalData[i,] <- sapply(strsplit(tempMatrix[i,],': '), function(x) x[[(length(x))]])
			}

			rownames(extraClinicalData) 					<- rownames(clinicalData)
			if(dataSet=='TCGA'){
				colnames(extraClinicalData) 				<- sapply(strsplit(tempMatrix[1,],': '), function(x) x[[1]])
			} else if(dataSet=='GSE30009'){
				colnames(extraClinicalData)[c(1:9,26:50)] 	<- sapply(strsplit(tempMatrix[1,c(1:9,26:50)],': '), function(x) x[[1]])
				colnames(extraClinicalData)[10:25] 			<- sapply(strsplit(tempMatrix[1,10:25],': '), function(x) x[[2]])
				colnames(extraClinicalData) 				<- gsub(' ','.',colnames(extraClinicalData))
				colnames(extraClinicalData) 				<- gsub('..','.',colnames(extraClinicalData),fixed=TRUE)
				colnames(extraClinicalData) 				<- gsub(';','',colnames(extraClinicalData))
				colnames(extraClinicalData) 				<- gsub('\\(','',colnames(extraClinicalData))
				colnames(extraClinicalData) 				<- gsub('\\)','',colnames(extraClinicalData))
				colnames(extraClinicalData)[which(colnames(extraClinicalData)=='surgical.debulking.or.residual.disease.cm')] <- 'residual_disease'
			} else if(dataSet=='GSE9891'){
				colnames(extraClinicalData)[c(1:9,15:55)] 	<- sapply(strsplit(tempMatrix[1,c(1:9,15:55)],': '), function(x) x[[1]])
				colnames(extraClinicalData)[10:14] 			<- sapply(strsplit(tempMatrix[1,10:14],': '), function(x) x[[2]])
			} else if(dataSet=='GSE32062.GPL6480'){
				colnames(extraClinicalData)[c(1:9,20:38)] 	<- sapply(strsplit(tempMatrix[1,c(1:9,20:38)],': '), function(x) x[[1]])
				colnames(extraClinicalData)[10:19] 			<- sapply(strsplit(tempMatrix[1,10:19],': '), function(x) x[[2]])
			} else if(dataSet=='GSE26712'){
				colnames(extraClinicalData)[c(1:9,14:37)] 	<- sapply(strsplit(tempMatrix[1,c(1:9,14:37)],': '), function(x) x[[1]])
				colnames(extraClinicalData)[10:13] 			<- sapply(strsplit(tempMatrix[1,10:13],': '), function(x) x[[2]])
			}

			if(any(extraClinicalData=='NA')) extraClinicalData[extraClinicalData=='NA'] <- NA

			renamedChosenExtraClinicalFeatures <- chosenExtraClinicalFeatures
			switch(dataSet,
				'GSE30009' 			= {	renamedChosenExtraClinicalFeatures <- gsub('tumor_residual_disease','residual_disease',renamedChosenExtraClinicalFeatures,fixed=TRUE)
										renamedChosenExtraClinicalFeatures <- gsub('preop_ca125','preop.ca125',renamedChosenExtraClinicalFeatures,fixed=TRUE)
										renamedChosenExtraClinicalFeatures <- gsub('platinum_sensitivity','platinum.sensitivity',renamedChosenExtraClinicalFeatures,fixed=TRUE)},
				'GSE9891' 			= {	renamedChosenExtraClinicalFeatures <- gsub('tumor_residual_disease','ResidualDisease',renamedChosenExtraClinicalFeatures,fixed=TRUE)
										renamedChosenExtraClinicalFeatures <- gsub('subtype','HistologicalSubtype',renamedChosenExtraClinicalFeatures,fixed=TRUE)
										renamedChosenExtraClinicalFeatures <- gsub('malignant_potential','Type',renamedChosenExtraClinicalFeatures,fixed=TRUE)},
				'GSE32062.GPL6480' 	= { renamedChosenExtraClinicalFeatures <- gsub("pfs (m)",'progression_free_survival_months',renamedChosenExtraClinicalFeatures,fixed=TRUE)
				        				renamedChosenExtraClinicalFeatures <- gsub("rec (1)",'recurrence',renamedChosenExtraClinicalFeatures,fixed=TRUE)
				        				renamedChosenExtraClinicalFeatures <- gsub("os (m)",'overall_survival_months',renamedChosenExtraClinicalFeatures,fixed=TRUE)
				        				renamedChosenExtraClinicalFeatures <- gsub("death (1)",'alive_dead',renamedChosenExtraClinicalFeatures,fixed=TRUE)}
			)
			extraClinicalFeatures <- as.data.frame(extraClinicalData,stringsAsFactors=FALSE)[,renamedChosenExtraClinicalFeatures,drop=FALSE]

			if(dataSet=='TCGA'&'tumor_residual_disease'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$tumor_residual_disease[which(extraClinicalFeatures$tumor_residual_disease=='No Macroscopic disease')]	<- 0
				extraClinicalFeatures$tumor_residual_disease[which(extraClinicalFeatures$tumor_residual_disease=='1-10 mm')] 				<- 1
				extraClinicalFeatures$tumor_residual_disease[which(extraClinicalFeatures$tumor_residual_disease=='11-20 mm')] 				<- 2
				extraClinicalFeatures$tumor_residual_disease[which(extraClinicalFeatures$tumor_residual_disease=='>20 mm')] 				<- 3
				extraClinicalFeatures$tumor_residual_disease 																				<- as.numeric(extraClinicalFeatures$tumor_residual_disease)
			}

			if(dataSet=='GSE30009'&'tumor_residual_disease'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$residual_disease[which(extraClinicalFeatures$residual_disease=='<1')]	<- 0
				extraClinicalFeatures$residual_disease[which(extraClinicalFeatures$residual_disease=='>1')] <- 1
				extraClinicalFeatures$residual_disease 														<- as.numeric(extraClinicalFeatures$residual_disease)
				extraClinicalFeatures$tumor_residual_disease 												<- extraClinicalFeatures$residual_disease
				extraClinicalFeatures$residual_disease 														<- NULL
			}

			if(dataSet=='GSE30009'&'chemotherapy'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$carboplatin 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$cisplatin 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$paclitaxel 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$gemcitabine 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$toxo 					<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$doxorubicin 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$irinotecan 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$tamoxifen 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$topotecan 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$cc2103 				<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$bevacizumab 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$capeciabine 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$anastrozole 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$cyclophosphamide 		<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$suboptTRSinterferon 	<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$letrozole 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$vinorelbine 			<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$carboplatin[grep("carboplatin|carbo|modifed triple doublets|tcA|platinum|t/c", extraClinicalFeatures$chemotherapy)] 	<- 1
				extraClinicalFeatures$cisplatin[grep("cisplatin", extraClinicalFeatures$chemotherapy)] 														<- 1
				extraClinicalFeatures$paclitaxel[grep("taxol|modifed triple doublets|tcA|taxane|t/c", extraClinicalFeatures$chemotherapy)] 					<- 1
				extraClinicalFeatures$gemcitabine[grep("gem|gemcitabine|Gemcitabine|modifed triple doublets", extraClinicalFeatures$chemotherapy)] 			<- 1
				extraClinicalFeatures$toxo[grep("toxo", extraClinicalFeatures$chemotherapy)] 																<- 1
				extraClinicalFeatures$doxorubicin[grep("doxil|modifed triple doublets|adriamyacin|doxorubicin", extraClinicalFeatures$chemotherapy)] 		<- 1
				extraClinicalFeatures$irinotecan[grep("cpt11", extraClinicalFeatures$chemotherapy)] 														<- 1
				extraClinicalFeatures$tamoxifen[grep("tamoxifen|tam", extraClinicalFeatures$chemotherapy)] 													<- 1
				extraClinicalFeatures$topotecan[grep("topo|topotecan|modifed triple doublets", extraClinicalFeatures$chemotherapy)] 						<- 1
				extraClinicalFeatures$cc2103[grep("cc2103", extraClinicalFeatures$chemotherapy)] 															<- 1
				extraClinicalFeatures$bevacizumab[grep("Avastin|tcA|avastin", extraClinicalFeatures$chemotherapy)] 											<- 1
				extraClinicalFeatures$capeciabine[grep("xeloda|Xelado", extraClinicalFeatures$chemotherapy)] 												<- 1
				extraClinicalFeatures$anastrozole[grep("arimidex", extraClinicalFeatures$chemotherapy)] 													<- 1
				extraClinicalFeatures$cyclophosphamide[grep("cytoxan", extraClinicalFeatures$chemotherapy)] 												<- 1
				extraClinicalFeatures$suboptTRSinterferon[grep("subopt tRS", extraClinicalFeatures$chemotherapy)] 											<- 1
				extraClinicalFeatures$letrozole[grep("letrozole|femara", extraClinicalFeatures$chemotherapy)] 												<- 1
				extraClinicalFeatures$vinorelbine[grep("navelbine", extraClinicalFeatures$chemotherapy)] 													<- 1
				extraClinicalFeatures$carboplatin 			<- as.numeric(extraClinicalFeatures$carboplatin)
				extraClinicalFeatures$cisplatin 			<- as.numeric(extraClinicalFeatures$cisplatin)
				extraClinicalFeatures$paclitaxel 			<- as.numeric(extraClinicalFeatures$paclitaxel)
				extraClinicalFeatures$gemcitabine 			<- as.numeric(extraClinicalFeatures$gemcitabine)
				extraClinicalFeatures$toxo 					<- as.numeric(extraClinicalFeatures$toxo)
				extraClinicalFeatures$doxorubicin 			<- as.numeric(extraClinicalFeatures$doxorubicin)
				extraClinicalFeatures$irinotecan 			<- as.numeric(extraClinicalFeatures$irinotecan)
				extraClinicalFeatures$tamoxifen 			<- as.numeric(extraClinicalFeatures$tamoxifen)
				extraClinicalFeatures$topotecan 			<- as.numeric(extraClinicalFeatures$topotecan)
				extraClinicalFeatures$cc2103 				<- as.numeric(extraClinicalFeatures$cc2103)
				extraClinicalFeatures$bevacizumab 			<- as.numeric(extraClinicalFeatures$bevacizumab)
				extraClinicalFeatures$capeciabine 			<- as.numeric(extraClinicalFeatures$capeciabine)
				extraClinicalFeatures$anastrozole 			<- as.numeric(extraClinicalFeatures$anastrozole)
				extraClinicalFeatures$cyclophosphamide 		<- as.numeric(extraClinicalFeatures$cyclophosphamide)
				extraClinicalFeatures$suboptTRSinterferon 	<- as.numeric(extraClinicalFeatures$suboptTRSinterferon)
				extraClinicalFeatures$letrozole 			<- as.numeric(extraClinicalFeatures$letrozole)
				extraClinicalFeatures$vinorelbine 			<- as.numeric(extraClinicalFeatures$vinorelbine)
				extraClinicalFeatures$chemotherapy 			<- NULL
			}

			if(dataSet=='GSE30009'&'response'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$response[which(extraClinicalFeatures$response=='PD')]			<- 0
				extraClinicalFeatures$response[which(extraClinicalFeatures$response=='SD')] 		<- 1
				extraClinicalFeatures$response[which(extraClinicalFeatures$response=='PR')] 		<- 2
				extraClinicalFeatures$response[which(extraClinicalFeatures$response=='CR')] 		<- 3
				extraClinicalFeatures$response[which(extraClinicalFeatures$response=='unknown')] 	<- NA
				extraClinicalFeatures$response 														<- as.numeric(extraClinicalFeatures$response)
			}

			if(dataSet=='GSE30009'&'platinum.sensitivity'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$platinum.sensitivity[which(extraClinicalFeatures$platinum.sensitivity=='sensitive')]	<- 1
				extraClinicalFeatures$platinum.sensitivity[which(extraClinicalFeatures$platinum.sensitivity=='resistant')] 	<- 0
				extraClinicalFeatures$platinum.sensitivity[which(extraClinicalFeatures$platinum.sensitivity=='unknown')] 	<- NA
				extraClinicalFeatures$platinum.sensitivity 																	<- as.numeric(extraClinicalFeatures$platinum.sensitivity)
			}

			if(dataSet=='GSE30009'&'preop.ca125'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$preop.ca125[which(extraClinicalFeatures$preop.ca125=='Elevated')]	<- NA
				extraClinicalFeatures$preop.ca125 <- as.numeric(extraClinicalFeatures$preop.ca125)
			}


			if(dataSet=='GSE9891'&'subtype'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$serous 																<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$endo 																	<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$adeno 																<- numeric(length(rownames(extraClinicalFeatures)))
				extraClinicalFeatures$serous[grep("Ser/PapSer", extraClinicalFeatures$HistologicalSubtype)] <- 1
				extraClinicalFeatures$endo[grep("Endo", extraClinicalFeatures$HistologicalSubtype)] 		<- 1
				extraClinicalFeatures$adeno[grep("Adeno", extraClinicalFeatures$HistologicalSubtype)] 		<- 1
				extraClinicalFeatures$HistologicalSubtype 													<- NULL
			}

			if(dataSet=='GSE9891'&'tumor_residual_disease'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$ResidualDisease[which(extraClinicalFeatures$ResidualDisease=='nil')]			<- 0
				extraClinicalFeatures$ResidualDisease[which(extraClinicalFeatures$ResidualDisease=='<1')]			<- 1
				extraClinicalFeatures$ResidualDisease[which(extraClinicalFeatures$ResidualDisease=='>1')] 			<- 2
				extraClinicalFeatures$ResidualDisease[which(extraClinicalFeatures$ResidualDisease=='NK')] 			<- NA
				extraClinicalFeatures$ResidualDisease[which(extraClinicalFeatures$ResidualDisease=='macrosizeNK')] 	<- NA
				extraClinicalFeatures$ResidualDisease 																<- as.numeric(extraClinicalFeatures$ResidualDisease)
				extraClinicalFeatures$tumor_residual_disease 														<- extraClinicalFeatures$ResidualDisease
				extraClinicalFeatures$ResidualDisease 																<- NULL
			}

			if(dataSet=='GSE9891'&'malignant_potential'%in%chosenExtraClinicalFeatures){
				extraClinicalFeatures$Type[which(extraClinicalFeatures$Type=='LMP')]	<- 0
				extraClinicalFeatures$Type[which(extraClinicalFeatures$Type=='MAL')]	<- 1
				extraClinicalFeatures$Type 												<- as.numeric(extraClinicalFeatures$Type)
				extraClinicalFeatures$malignant_potential 								<- extraClinicalFeatures$Type
				extraClinicalFeatures$Type 												<- NULL
			}
		} else extraClinicalFeatures <- NA
	} else {
		clinicalFeatures <- clinicalData
		names(clinicalFeatures)[which(names(clinicalFeatures)=='vital_status')] 						<- 'event'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='days_to_death')] 						<- 'days_to_event'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='age_at_initial_pathologic_diagnosis')] 	<- 'age'
		names(clinicalFeatures)[which(names(clinicalFeatures)=='days_to_tumor_recurrence')] 			<- 'days_to_recurrence'
		clinicalFeatures <- clinicalFeatures[,c('event','days_to_event')]
	}

	if((clinicalFeaturesFlag[1]&clinicalFeaturesFlag[2]&geneExpressionFlag)|(!clinicalFeaturesFlag[1]&clinicalFeaturesFlag[2]&geneExpressionFlag)){
		inputData <- merge(clinicalFeatures,extraClinicalFeatures,by='row.names')
		inputData <- merge(inputData,expressionFeatures,by.x='Row.names',by.y='row.names')
	} else if((clinicalFeaturesFlag[1]&!clinicalFeaturesFlag[2]&geneExpressionFlag)|(!clinicalFeaturesFlag[1]&!clinicalFeaturesFlag[2]&geneExpressionFlag)){
		inputData <- merge(clinicalFeatures,expressionFeatures,by='row.names')
	} else if((clinicalFeaturesFlag[1]&clinicalFeaturesFlag[2]&!geneExpressionFlag)|(!clinicalFeaturesFlag[1]&!clinicalFeaturesFlag[2]&!geneExpressionFlag)){
		inputData <- merge(clinicalFeatures,extraClinicalFeatures,by='row.names')
	} else if (clinicalFeaturesFlag[1]&!clinicalFeaturesFlag[2]&!geneExpressionFlag){
		inputData <- clinicalFeatures
	}

	if(length(which(is.na(inputData$'days_to_event')))!=0) inputData = inputData[-which(is.na(inputData$'days_to_event')),]
	inputData$event = sapply(1:dim(inputData)[1], function(x) if(inputData$event[x]=='deceased'){1} else {0})
	names(inputData)[which(names(inputData)=='Row.names')] = 'patient_num'
	rownames(inputData) = inputData$patient_num

	chosenModelFeatures = names(inputData)[-which(names(inputData) %in% c('patient_num','days_to_recurrence','days_to_event','event'))]

	toReturn = list('inputData'=inputData,'allFeatureNames'=chosenModelFeatures)
	if(clinicalFeaturesFlag[1]|clinicalFeaturesFlag[2]){
		toReturn$clinicalFeatures 			<- clinicalFeatures
		toReturn$extraClinicalFeatures 		<- extraClinicalFeatures
		toReturn$clinicalFeatureNames 		<- chosenClinicalFeatures
		toReturn$extraClinicalFeatureNames 	<- chosenExtraClinicalFeatures
	}
	if(geneExpressionFlag){
		toReturn$expressionFeatures 		<- expressionFeatures
		toReturn$expressionFeatureNames 	<- chosenExpressionFeatures
	}

return(toReturn)

}