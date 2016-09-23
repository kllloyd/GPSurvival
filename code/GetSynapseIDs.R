GetSynapseIDs <- function(parameterStructure){
	#-------------------------------------------------------------------------------------------------------#
    # K Lloyd 2016_09_16
    #-------------------------------------------------------------------------------------------------------#

	cancer 			<- parameterStructure$cancer  	# c('KIRC','OV','GBM','LUSC')
	molPlatform 	<- parameterStructure$molPlatform 	# c('SCNA','methyl','mRNA','miRNA','protein','None')
	clinicalFlag 	<- parameterStructure$clinicalFlag

	switch(cancer,
		'KIRC' 	= {switch(molPlatform,
					'SCNA'		={exp.ID <- 'syn1710287'},
					'methyl'	={exp.ID <- 'syn1710289'},
					'mRNA'		={exp.ID <- 'syn1710293'},
					'miRNA'		={exp.ID <- 'syn1710291'},
					'protein'	={exp.ID <- 'syn1710306'},
					'None'		={exp.ID <- ''}
					)
					if(clinicalFlag) clinical.ID <- 'syn1715824' else clinical.ID <- ''
					surv.ID 	<- 'syn1710303'
					train.ID 	<- 'syn1714093'
					test.ID 	<- 'syn1714090'
				},
		'OV' 	= {switch(molPlatform,
					'SCNA'		={exp.ID <- 'syn1710316'},
					'methyl'	={exp.ID <- 'syn1710320'},
					'mRNA'		={exp.ID <- 'syn1710361'},
					'miRNA'		={exp.ID <- 'syn1710359'},
					'protein'	={exp.ID <- 'syn1710314'},
					'None'		={exp.ID <- ''}
					)
					if(clinicalFlag) clinical.ID <- 'syn1715828' else clinical.ID <- ''
					surv.ID 	<- 'syn1710363'
					train.ID 	<- 'syn1714105'
					test.ID 	<- 'syn1714102'
				},
		'GBM' 	= {switch(molPlatform,
					'SCNA'		={exp.ID <- 'syn1710366'},
					'methyl'	={exp.ID <- 'syn1710374'},
					'mRNA'		={exp.ID <- 'syn1710372'},
					'miRNA'		={exp.ID <- 'syn1710368'},
					'protein'	={exp.ID <- ''},
					'None'		={exp.ID <- ''}
					)
					if(clinicalFlag) clinical.ID <- 'syn1715822' else clinical.ID <- ''
					surv.ID 	<- 'syn1710370'
					train.ID 	<- 'syn1714087'
					test.ID 	<- 'syn1714083'
				},
		'LUSC' 	= {switch(molPlatform,
					'SCNA'		={exp.ID <- 'syn1710378'},
					'methyl'	={exp.ID <- ''},
					'mRNA'		={exp.ID <- 'syn1710382'},
					'miRNA'		={exp.ID <- 'syn1710380'},
					'protein'	={exp.ID <- 'syn1710386'},
					'None'		={exp.ID <- ''}
					)
					if(clinicalFlag) clinical.ID <- 'syn1715826' else clinical.ID <- ''
					surv.ID 	<- 'syn1710384'
					train.ID 	<- 'syn1714099'
					test.ID 	<- 'syn1714096'
				}
		)

	toReturn 				<- parameterStructure
	toReturn$exp.ID 		<- exp.ID
	toReturn$clinical.ID 	<- clinical.ID
	toReturn$surv.ID 		<- surv.ID
	toReturn$train.ID 		<- train.ID
	toReturn$test.ID 		<- test.ID

	return(toReturn)
}