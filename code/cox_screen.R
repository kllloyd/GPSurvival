# Yuan Yuan 2012-12-13
# This script calculates screen for the significant features against the clinical outcome, using univariate cox model
library(foreach)
# library(doMC)
# registerDoMC(6)

cox.screen <- function(y.train,x.train, pvalue.cutoff=0.05, qvalue.cutoff=0.2, top=100) 
  {
    # calculate the p-value from likelihood ratio test of univariate cox model
    # for categorical variable, such as mutation, should use log-rank test
    feature.pvalue <- c()
    feature.name <- c()
    feature.col <- c()
        #for (j in 1:ncol(x.train))
    result <- foreach (j= 1:ncol(x.train),.errorhandling='remove',.combine=rbind) %dopar%
    {
      feature <- colnames(x.train)[j]
      x <- x.train[,j]
      if (length(which(table(x)> 0.8*length(x)))>0) # discard the flat values, e.g. zeros for RNAseq and miRNAseq
        {
          stop()
        }
      cox <- try(coxph(y.train~x))
      if (class(cox)=="try-error")
        {
          #print(j)
          stop() #next# note next does not work with foreach %dopar%    
        }
      p.value <-summary(cox)$logtest["pvalue"] # likelihood ratio test
      list(p.value, feature, j)
    }
       
    feature.pvalue <- unlist(result[,1])
    print(paste("Total number of valid features (after removal of potential flat records):", length(feature.pvalue)))
    feature.name <- unlist(result[,2])
    feature.col <- unlist(result[,3])
          
    names(feature.pvalue)=c()
    names(feature.name)=c()
    names(feature.col)=c()
    
    q.value <- p.adjust(feature.pvalue, method="fdr")
    for (i in 1:10/10)
      {
        print(paste("qvalue <=", i, ":", length(which(q.value<=i))))
      }
    col.sig <- which( feature.pvalue < pvalue.cutoff) # & q.value < qvalue.cutoff  )
    print(paste("Significant records: p-value < 0.05:", length(col.sig)))  #and FDR <",qvalue.cutoff,":", length(col.sig)))
   
    name.sig <- feature.name[col.sig]
    pvalue.sig <- feature.pvalue[col.sig]
    #qvalue.sig <- q.value[col.sig]
    

    col.retain <- feature.col[col.sig]
    if (length(col.retain)> top) # only keep the top significant ones if there are too many
      {
        col.retain <- col.retain[head(sort(pvalue.sig, index.return=T)$ix, n=top)]
      }
    
    return(col.retain)
}
