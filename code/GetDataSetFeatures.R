GetDataSetFeatures <- function(dataSource,dataSet,geneSubsetFlag,clinicalSubsetFlag){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#
	# geneSubsetFlag: 		'All' 		all genes in data set
	#						'SRGS1' 	genes from systematic review selected by at least 2 studies
	# 						'SRGS2'		genes from systematic review selected by at least 3 studies
	#						'TaqMan'	genes included on OvCa TaqMan array
	# 						'MP3'		genes idenfified as informative from TaqMan by Masters MiniProject 3
	#						'None' 		no genes
	# Clinical features are highly dependant on the data set in use. 
	#-------------------------------------------------------------------------------------------------------#

	switch(geneSubsetFlag,
		'All'	={chosenExpressionFeatures 	<- 'All'},
		'SRGS1'	={chosenExpressionFeatures 	<- c('AGR2','MUTYH','AKAP12','TP53','TOP2A','FOXA2','SRC','SIVA1','ALDH9A1','LGR5','EHF','BAX','CES2','CPE','FGFBP1','TUBB4A',
												'ZNF12','RBM39','RFC3','GNPDA1','ANXA3','NFIB','ACTR3B','YWHAE','CYP51A1','HMGCS1','ZMYND11','FADS2','SNX7','ARHGDIA',
												'NDST1','DAP','ERCC8','GUCY1B3','HDAC1','HDAC2','IGFBP5','IL6','LSAMP','DGKZ','MYCBP','S100A10','SLC1A3','NCOA1','TIAM1',
												'VEGFA','RPL36','LBR','ABCB1','FASLG','TIMP1','FN1','TGFB1','XPA','POLH','ITGAE','ZNF200','COL3A1','CXCR7','EPHB3','NBN',
												'PCF11','DFNB31','BRCA2','AADAC','CD38','CHIT1','CXCR4','EFNB2','MECOM','FILIP1L','HSPB7','LRIG1','MMP1','PSAT1','SDF2L1',
												'TCF15','EPHB2','ETS1','TRIM27','MARK4','B4GALT5','ABCB10','AOC1')},
		'SRGS2'	={chosenExpressionFeatures 	<- c('MUTYH','AKAP12','TP53','TOP2A','FOXA2','AGR2')},
		'TaqMan'={chosenExpressionFeatures 	<- c('AKT1','AKT2','AKT3','APAF1','BAD','BAX','BCL2','BCL2L1','BID','CFLAR','FAS','FASLG','HSPD1','HSPA1A','HSPA1L','HSP90AA1',
												'HSP90AB1','HSP90B1','BIRC2','IGF1','IGF1R','IGF2','IGF2R','IGFBP1','IGFBP2','DNAJC15','MCL1','FRAP1','NFKB1','PIK3CA',
												'PTEN','STAT3','BIRC5','BIRC4','ATP7B','ABCG2','CES1','CES2','NT5C2','DPYD','FPGS','H2AFX','GCLC','GCLM','GSTP1','SLC29A1',
												'SLC29A2','ABCB1','ABCC1','ABCC2','ABCC3','ABCC4','ABCC5','ABCC6','ABCC8','MVP','UMPS','RRM1','SOD1','TAP1','TAP2','ABCB4',
												'TYMS','HPRT1','HMBS','SDHA','TBP','ATM','BRCA1','ERCC1','ERCC2','MGMT','MLH1','MSH2','MSH6','RAD51','TOP1','TOP2A','TOP2B',
												'XPA','XRCC1','XRCC5','XRCC6','APC','TUBB3','PTGS2','EGFR','ERBB2','ERBB3','ERBB4','HIF1A','MKI67','CDKN2A','CDKN1A',
												'CDKN1B','TP53','VEGFA')},
		'MP3'	={chosenExpressionFeatures 	<- c('ABCC2','BIRC5','CA2','CASP1','DCK','ERAL1','ERBB3','ESR2','GADD45B','GSTM1','PGR','POLE','PTEN','RAD51','SDHA',
												'SLC28A3','VEGFA')},
		'None'	={chosenExpressionFeatures 	<- 'None'}
	)

	if(dataSource=='CuratedOvarian'){
		switch(dataSet,
			'GSE30009' 			={switch(clinicalSubsetFlag,
									'None' 	={chosenClinicalFeatures 		<- 'None'
											  chosenExtraClinicalFeatures 	<- 'None'},
									'One'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','debulking','pltx','tax',
						  														'days_to_death','vital_status')
						  					  chosenExtraClinicalFeatures 	<- 'None'},
									'Two'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','debulking','pltx','tax',
						  														'days_to_death','vital_status')
						  					  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','preop_ca125','platinum_sensitivity')},
									'Three'	={chosenClinicalFeatures 		<- 'None'
						  				  	  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','preop_ca125','platinum_sensitivity')},
						  			'Four'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death',
						  														'vital_status')
						  					  chosenExtraClinicalFeatures 	<- 'None'},
						  			'Five'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death','debulking',
						  														'vital_status')
						  					  chosenExtraClinicalFeatures 	<- 'None'},
									'Five'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death','debulking',
						  														'vital_status')
						  					  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','preop_ca125')})},
			'GSE9891'			={switch(clinicalSubsetFlag,
								 	'None' 	={chosenClinicalFeatures 		<- 'None'
											  chosenExtraClinicalFeatures 	<- 'None'},
									'One'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','debulking','pltx','tax',
								 												'days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
									'Two'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','debulking','pltx','tax',
																				'days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','subtype','malignant_potential')},
									'Three'	={chosenClinicalFeatures 		<- 'None'
								 		 	  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','subtype','malignant_potential')},
								  	'Four'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death',
								  												'vital_status')
								  			  chosenExtraClinicalFeatures 	<- 'None'},
								  	'Five'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death','debulking',
								  												'vital_status')
								  			  chosenExtraClinicalFeatures 	<- 'None'},
								  	'Six'	={chosenClinicalFeatures 		<- c('grade','tumorstage','age_at_initial_pathologic_diagnosis','days_to_death',
								  												'vital_status')
								  			  chosenExtraClinicalFeatures 	<- c('tumor_residual_disease','malignant_potential')})},
			'GSE32062.GPL6480' 	={switch(clinicalSubsetFlag, 
									'None'	={chosenClinicalFeatures 		<- 'None'
											  chosenExtraClinicalFeatures 	<- 'None'},
									'One'	={chosenClinicalFeatures 		<- c('grade','tumorstage','debulking','pltx','tax','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
									'Two'	={chosenClinicalFeatures 		<- c('grade','tumorstage','debulking','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
								 	'Three'	={chosenClinicalFeatures 		<- c('grade','tumorstage','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'})},
			'GSE26712' 			={switch(clinicalSubsetFlag, 
									'None'	={chosenClinicalFeatures 		<- 'None'
											  chosenExtraClinicalFeatures 	<- 'None'},
									'One'	={chosenClinicalFeatures 		<- c('tumorstage','recurrence_status','debulking','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
								 	'Two'	={chosenClinicalFeatures 		<- c('tumorstage','substage','debulking','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
									'Three'	={chosenClinicalFeatures 		<- c('tumorstage','debulking','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'},
								 	'Four'	={chosenClinicalFeatures 		<- c('tumorstage','substage','days_to_death','vital_status')
								 			  chosenExtraClinicalFeatures 	<- 'None'})}
		)
	} else if(dataSource=='TCGA2STAT'){
		switch(clinicalSubsetFlag,
			'None' 	={chosenClinicalFeatures <- 'None'},
			'One'	={chosenClinicalFeatures <- 'yearstobirth'},
			'Two'	={chosenClinicalFeatures <- 'pathologicalstage'},
			'Three' ={chosenClinicalFeatures <- 'residualtumor'},
			'Four'	={chosenClinicalFeatures <- c('yearstobirth','pathologicalstage')},
			'Five'	={chosenClinicalFeatures <- c('yearstobirth','pathologicalstage','residualtumor')},
			'Six' 	={chosenClinicalFeatures <- c('pathologicalstage','residualtumor')}
		)
	} else if(dataSource=='TCGASynapse'){
		switch(dataSet[1],
			'BRCA' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('age_at_initial_pathologic_diagnosis',
															  'year_of_initial_pathologic_diagnosis',
															  'ajcc_neoplasm_disease_stage',
															  'breast_carcinoma_progesterone_receptor_status',
															  'breast_carcinoma_estrogen_receptor_status')},
						'Two'	={chosenClinicalFeatures <- c('age_at_initial_pathologic_diagnosis',
															  'year_of_initial_pathologic_diagnosis')},
						'Three' ={chosenClinicalFeatures <- c('ajcc_neoplasm_disease_stage',
															  'breast_carcinoma_progesterone_receptor_status',
															  'breast_carcinoma_estrogen_receptor_status')})},
			'LUAD' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('gender',
															  'tumor_stage',
															  'year_of_initial_pathologic_diagnosis',
															  'days_to_birth',
															  'tobacco_smoking_history_indicator')},
						'Two'	={chosenClinicalFeatures <- c('gender',
															  'tumor_stage',
															  'days_to_birth',
															  'tobacco_smoking_history_indicator')},
						'Three' ={chosenClinicalFeatures <- c('tumor_stage',
															  'days_to_birth',
															  'tobacco_smoking_history_indicator')})},
			'THCA' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('year_of_initial_pathologic_diagnosis',
															  'ajcc_neoplasm_disease_stage',
															  'days_to_birth',
															  'gender',
															  'residual_tumor')},
						'Two'	={chosenClinicalFeatures <- c('year_of_initial_pathologic_diagnosis',
															  'ajcc_neoplasm_disease_stage',
															  'days_to_birth',
															  'gender')},
						'Three' ={chosenClinicalFeatures <- c('ajcc_neoplasm_disease_stage',
															  'days_to_birth',
															  'gender')})},
			'HNSC' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('days_to_birth',
															  'neoplasm_histologic_grade',
															  'gender',
															  'tobacco_smoking_history_indicator',
															  'year_of_initial_pathologic_diagnosis')},
						'Two'	={chosenClinicalFeatures <- c('days_to_birth',
															  'neoplasm_histologic_grade',
															  'gender',
															  'tobacco_smoking_history_indicator')},
						'Three' ={chosenClinicalFeatures <- c('days_to_birth',
															  'neoplasm_histologic_grade',
															  'tobacco_smoking_history_indicator')})},
			'GBM' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('gender',
															  'prior_glioma',
															  'age_at_initial_pathologic_diagnosis',
															  'year_of_initial_pathologic_diagnosis',
															  'karnofsky_performance_score')},
						'Two'	={chosenClinicalFeatures <- c('prior_glioma',
															  'age_at_initial_pathologic_diagnosis',
															  'karnofsky_performance_score')},
						'Three' ={chosenClinicalFeatures <- c('prior_glioma',
															  'age_at_initial_pathologic_diagnosis')})},
			'OV' 	={switch(clinicalSubsetFlag,
						'None' 	={chosenClinicalFeatures <- 'None'},
						'One'	={chosenClinicalFeatures <- c('neoplasm_histologic_grade',
															  'tumor_residual_disease',
															  'tumor_stage',
															  'age_at_initial_pathologic_diagnosis',
															  'year_of_initial_pathologic_diagnosis')},
						'Two'	={chosenClinicalFeatures <- c('neoplasm_histologic_grade',
															  'tumor_stage',
															  'age_at_initial_pathologic_diagnosis',
															  'year_of_initial_pathologic_diagnosis')},
						'Three' ={chosenClinicalFeatures <- c('neoplasm_histologic_grade',
															  'tumor_stage',
															  'age_at_initial_pathologic_diagnosis')})}
		)
	}
	
	toReturn <- list('chosenExpressionFeatures'=chosenExpressionFeatures,'chosenClinicalFeatures'=chosenClinicalFeatures)
	if(dataSource=='CuratedOvarian') toReturn$chosenExtraClinicalFeatures <- chosenExtraClinicalFeatures
	return(toReturn)
}