#'Sample-specific differentially expressed (DE) genes
#'
#'This function computes the median expression of a gene across all samples and further computed the fold change of that gene to the computed median expression. A gene was considered up if the fold-change ≥ 4 and down if the fold-change ≤ −4.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 3 | Function 1 |
#'
#' @param RnaSeq_data A R dataframe with gene-centric expression values of bulk tumor samples. Rows are the genes, columns are the tumor samples.
#' @param parallel A boolean value ('TRUE' or 'FALSE')indicating if the users want to run this function in multi-core.
#' @param ncores An integer value for the number of cores.
#' @param CTDDirectory A character string for the file path of the directory for the output to be stored.
#'
#' @return R dataframe objects each with a single column containing DE genes for each tumor sample. The column is labeled as gene.
#' @export
#'
#' @examples
GetSampleDEFun<-function(RnaSeq_data,parallel=c('TRUE','FALSE'),ncores=2, CTDDirectory="~"){
  library(matrixStats)
  library(parallel)
  library(pbapply)
  if((is(RnaSeq_data, "SummarizedExperiment"))){
    RnaSeq_data<-assay(RnaSeq_data)
  }

  dir.create(paste0(CTDDirectory,"/Patient_DE_genes/"))
  dir.create(paste0(CTDDirectory,"/Patient_DE_genes/UP_Gene/"))
  dir.create(paste0(CTDDirectory,"/Patient_DE_genes/DOWN_Gene/"))
  dir.create(paste0(CTDDirectory,"/Patient_DE_genes/DE_Gene/"))
  processFiles <- function(fileNames){

    SampleName <- fileNames
    RnaSeq_data<-as.data.frame(RnaSeq_data)
    gene_expr<-RnaSeq_data[, SampleName,drop=FALSE]
    gene_med<-RnaSeq_data[,c(1,2),drop=FALSE]
    gene_expr_med<-cbind(gene_med,gene_expr)
    colnames(gene_expr_med)<-c("gene","median_expr","expr")
    gene_expr_med[gene_expr_med$median_expr==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
    gene_expr_med[gene_expr_med$expr==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
    gene_expr_med$fold_change<-gene_expr_med$expr/gene_expr_med$median_expr
    upgene<-gene_expr_med[gene_expr_med$fold_change>=4,]
    downgene<-gene_expr_med[gene_expr_med$fold_change<=0.25,]
    pat_de_genes <- data.table(rbind(upgene,downgene))[,1]

    save(upgene,file=paste0(CTDDirectory,"/Patient_DE_genes/UP_Gene/",SampleName,"_upgene.Rda"))
    save(downgene,file=paste0(CTDDirectory,"/Patient_DE_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
    save(pat_de_genes,file=paste0(CTDDirectory,"/Patient_DE_genes/DE_Gene/",SampleName,"_pat_de_genes.Rda"))

  }


  x <- data.table(rowMedians(RnaSeq_data,na.rm=TRUE))
  colnames(x)<-"Median_expr"
  row.names(x)<-c(rownames(RnaSeq_data))
  Rnames<-rownames(RnaSeq_data)
  RnaSeq_data<-cbind(x,RnaSeq_data)
  RnaSeq_data<-as.data.frame(RnaSeq_data)
  RnaSeq_data<-cbind(Rnames, RnaSeq_data)

  n<-ncol(RnaSeq_data)
  n<-as.integer(n)
  fileNames1<-colnames(RnaSeq_data)
  fileNames<-fileNames1[3:n]

  if (parallel==TRUE){
    print("cluster making started")
    cl <- makeCluster(mc <- getOption("cl.cores", ncores),type = "PSOCK")
    invisible(clusterEvalQ(cl, library('glmnet')))
    invisible(clusterEvalQ(cl, library(janitor)))
    invisible(clusterEvalQ(cl, library(stringr)))
    invisible(clusterEvalQ(cl, library(plyr)))
    invisible(clusterEvalQ(cl, library(data.table)))
    invisible(clusterEvalQ(cl, library(readr)))
    invisible(clusterEvalQ(cl, library(gsubfn)))
    invisible(clusterEvalQ(cl, library(parallel)))
    invisible(clusterEvalQ(cl, library(sqldf)))
    invisible(clusterEvalQ(cl, library(reshape)))
    invisible(clusterEvalQ(cl, library(reshape2)))
    invisible(clusterEvalQ(cl, library(pbapply)))

    print("cluster export started")
    parallel::clusterExport(cl=cl,c("RnaSeq_data"),envir=environment())
    print("cluster making done")

    parLapplyLB(cl,fileNames,processFiles)

    stopCluster(cl)
    print("cluster stopped")
  }
  if (parallel==FALSE){
    lapply(fileNames,processFiles)
  }
  return()
}
