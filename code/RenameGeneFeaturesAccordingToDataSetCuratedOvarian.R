RenameGeneFeaturesAccordingToDataSetCuratedOvarian <- function(dataSetExpression,chosenExpressionFeatures){
	#--------------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #--------------------------------------------------------------------------------------------------------------#
    # List of gene synonyms "protein-coding_gene.txt" downloaded from http://www.genenames.org/cgi-bin/statistics
    #--------------------------------------------------------------------------------------------------------------#

	## Set up list of gene synonyms from HGNC database ##
	hgncDatabase 		<- read.delim("protein-coding_gene.txt",stringsAsFactors=FALSE)
	geneListHGNC 		<- vector("list", dim(hgncDatabase)[1])
	names(geneListHGNC) <- gsub("[[:punct:]]", "", hgncDatabase$Approved.Symbol)
	for (i in 1:length(geneListHGNC)){
	  geneListHGNC[[names(geneListHGNC)[i]]] <- c(hgncDatabase$Approved.Symbol[i],strsplit(hgncDatabase$Previous.Symbols[i], ', ', fixed = TRUE)[[1]],strsplit(hgncDatabase$Synonyms[i], ', ', fixed = TRUE)[[1]])
	  geneListHGNC[[names(geneListHGNC)[i]]] <- geneListHGNC[[names(geneListHGNC)[i]]][geneListHGNC[[names(geneListHGNC)[i]]]!= ""]
	}

	## Set up list of genes in data source ##
	geneListDataSource 		<- row.names(dataSetExpression)
	geneListDataSource 		<- lapply(geneListDataSource,function(ch) ifelse(length(grep('///',ch,fixed=TRUE))!=0,strsplit(ch,'///',fixed=TRUE),ch))

	# Set up list of genes wanted as features ##
	geneListChosenFeatures 	<- chosenExpressionFeatures
	
	## Find genes in geneListChosenFeatures present in geneListDataSet ##
	geneListRenamed 		<- geneListChosenFeatures
	names(geneListRenamed) 	<- geneListChosenFeatures
	for(i in 1:length(geneListChosenFeatures)){
		geneListCompare 	<- sapply(1:length(geneListDataSource),function(x) c(ifelse(any(geneListChosenFeatures[i]%in%geneListDataSource[[x]][[1]]),geneListChosenFeatures[i],NA),ifelse(any(geneListChosenFeatures[i]%in%geneListDataSource[[x]][[1]]),paste0(geneListDataSource[[x]][[1]],collapse='///'),NA)))
		geneListRenamed[i] 	<- ifelse(length(geneListCompare[2,which(geneListCompare[1,]==geneListRenamed[i])])!=0,geneListCompare[2,which(geneListCompare[1,]==geneListRenamed[i])],NA)
	}

	## Find genes missing from geneListDataSet ##
	geneListMissing <- names(geneListRenamed[is.na(geneListRenamed)])
	
	if(length(geneListMissing)!=0){
		## Find alternative names for missing genes in geneListHGNC ##
		geneListAltNames <- list()
		for(i in 1:length(geneListMissing)){
			geneListAltNames[[i]] <- sapply(1:length(geneListHGNC),function(x) if(any(geneListMissing[i]%in%geneListHGNC[[x]])) geneListHGNC[[x]] else NA)
			geneListAltNames[[i]] <- geneListAltNames[[i]][!is.na(geneListAltNames[[i]])]
			names(geneListAltNames)[[i]] <- geneListMissing[i]
		}
		## Find genes in list of alternative names in geneListDataSet ##
		geneListFound <- character()
		for(i in 1:length(geneListAltNames)){
			if(is.list(geneListAltNames[[i]])){
				geneListCompare <- sapply(1:length(geneListDataSource),function(y) ifelse(any(geneListAltNames[[i]][[1]]%in%geneListDataSource[[y]][[1]]),paste0(geneListDataSource[y],collapse='///'),NA))
			} else {
				geneListCompare <- rep(NA,length(geneListDataSource))
			}
			geneListFound[i] <- ifelse(any(!is.na(geneListCompare)),geneListCompare[!is.na(geneListCompare)],NA)
			names(geneListFound)[i] <- geneListMissing[i]
		}
		geneListRenamed[names(geneListRenamed)%in%names(geneListFound)] <- geneListFound 	# Any genes that have not been found are marked as NA
	}

	## Tidy up and prepare to return ##
	geneListPresentInDataSet <- geneListRenamed[!is.na(geneListRenamed)]				# Version of geneListRenamed without NA genes
	
	nGenesInDataSet <- sum(!is.na(geneListRenamed))
	
	toReturn <- list('geneListRenamed'=geneListRenamed,'geneListPresentInDataSet'=geneListPresentInDataSet,'nGenesInDataSet'=nGenesInDataSet)
return(toReturn)
}
