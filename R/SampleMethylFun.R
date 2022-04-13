#'Sample-specific deconvoluted DNA methylation profile.
#'
#'This function utilizes the estimated proportion of constituent cell types and the estimated gene-centric DNA methylation profile of constituent cell type to compute the sample-specific deconvoluted DNA methylation profile for each tumor sample.
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 1 | Function 4 |
#' @param Deconv_meth_gene A R dataframe with the deconvoluted DNA methylation values for each gene in constituent cell types.
#' @param Deconv_proportions A R dataframe with the estimated proportions of constituent cell types in samples.
#' @param CTDDirectory A character string for the file path of the directory for the output to be stored.
#'
#' @return R dataframes including the deconvoluted DNA methylation profile for each sample in the specified directory. The rows are the genes and the columns are the cell types.This function labels the columns as V1, V2, V3,....etc.
#' @export
#'
#' @examples
SampleMethylFun<-function(Deconv_meth_gene,Deconv_proportions, CTDDirectory="~"){

  if((is(Deconv_meth_gene, "SummarizedExperiment"))){
    Deconv_meth_gene<-assay(Deconv_meth_gene)
  }

  if((is(Deconv_proportions, "SummarizedExperiment"))){
    Deconv_proportions<-assay(Deconv_proportions)
  }

  Deconv_meth_gene<-abs(Deconv_meth_gene)
  Deconv_proportions<-abs(Deconv_proportions)
  dir.create(paste0(CTDDirectory, "/Sample_methylation"))

  processMethExp<- function(sample){
    tryCatch(
      {
        #sample=samples$V1[5]
        proportions<-Deconv_proportions[unlist(sample),]
        n1<-ncol(Deconv_meth_gene)
        if(n1==3){
          methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                             Deconv_meth_gene[,3]*proportions[3])
        }
        if(n1==4){
          methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                             Deconv_meth_gene[,3]*proportions[3],Deconv_meth_gene[,4]*proportions[4])
        }
        if(n1==5){
          methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                             Deconv_meth_gene[,3]*proportions[3],Deconv_meth_gene[,4]*proportions[4],
                             Deconv_meth_gene[,5]*proportions[5])
        }

        methylation<-as.matrix(methylation)
        colnames(methylation)<-colnames(Deconv_meth_gene)
        row.names(methylation)<-row.names(Deconv_meth_gene)
        save(methylation,file=paste0(CTDDirectory,"/Sample_methylation/",sample,".Rda"))


      }, error = function(error_condition) {
        cat(paste0(sample," : ",error_condition),file="~/processMethExp.txt",sep="\n",append=TRUE)
        return()
      }, finally={
      }
    )}

  samples<-data.table(row.names(Deconv_proportions))
  lapply(samples$V1,processMethExp)
  return()
}

