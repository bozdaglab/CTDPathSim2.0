#'Cell line-specific differentially (DM) genes
#'
#'This function computes cell line-specific DM genes for each cell lines It computes M-values from gene-centric beta values and computes the median of gene-centric M-values across all the cell lines It considers a gene hypermethylated if the M-value fold-change ≥ 4 and hypomethylated if fold-change ≤ -4.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 4 | Function 2 |
#'
#' @param DNAmethylation_data A R dataframe with gene-centric DNA methylation beta values of cell lines. Rows are the genes, columns are the cell lines.
#' @param parallel A boolean value ('TRUE' or 'FALSE')indicating if the users want to run this function in multi-core.
#' @param ncores A integer value for the number of cores.
#' @param CTDDirectory A character string for the file path of the directory for the output to be stored.
#'
#' @return R dataframe objects each with a single column containing DM genes for each cell line. The column is labeled as gene.
#' @export
#'
#' @examples
GetCellLineDMFun<-function(DNAmethylation_data,parallel=c('TRUE','FALSE'),ncores=2, CTDDirectory="~"){
  if((is(DNAmethylation_data, "SummarizedExperiment"))){
    DNAmethylation_data<-assay(DNAmethylation_data)
  }
  dir.create(paste0(CTDDirectory,"/CellLine_DM_genes/"))
  dir.create(paste0(CTDDirectory,"/CellLine_DM_genes/UP_Gene/"))
  dir.create(paste0(CTDDirectory,"/CellLine_DM_genes/DOWN_Gene/"))
  dir.create(paste0(CTDDirectory,"/CellLine_DM_genes/DM_Gene/"))
  processFiles <- function(fileNames){
    SampleName <- fileNames
    DNAmethylation_data<-as.data.frame(DNAmethylation_data)
    gene_meth<-DNAmethylation_data[,SampleName,drop=FALSE]
    gene_med<-DNAmethylation_data[,c(1,n),drop=FALSE]
    gene_meth_med<-cbind(gene_med,gene_meth)
    colnames(gene_meth_med)<-c("gene","median_meth","meth")
    gene_meth_med[gene_meth_med$median_meth==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
    gene_meth_med[gene_meth_med$meth==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
    gene_meth_med$fold_change<-abs(gene_meth_med$meth/gene_meth_med$median_meth)
    upgene<-gene_meth_med[gene_meth_med$fold_change>=4,]
    downgene<-gene_meth_med[gene_meth_med$fold_change<=0.25,]
    pat_dm_genes <- data.table(rbind(upgene,downgene))[,1]


    save(upgene,file=paste0(CTDDirectory,"/CellLine_DM_genes/UP_Gene/",SampleName,"_upgene.Rda"))
    save(downgene,file=paste0(CTDDirectory,"/CellLine_DM_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
    save(pat_dm_genes,file=paste0(CTDDirectory,"/CellLine_DM_genes/DM_Gene/",SampleName,"_cell_dm_genes.Rda"))

  }


  DNAmethylation_data<-as.data.frame(DNAmethylation_data)
  DNAmethylation_data[DNAmethylation_data==0]<-0.00000000001 ## to avoid the log transformation of 0
  #Writing a function for transforming beta values to M values
  Mtranform <- function(x) {log(x[1]/(1-x[1]))}
  n<-ncol(DNAmethylation_data)
  DNAmethylation_data[, seq_len(n) ] <- as.data.frame(lapply(DNAmethylation_data[, seq_len(n) ], FUN = function(x) {vapply(x, FUN = Mtranform, FUN.VALUE = 1.0)}))
  x <- data.table(rowMedians(as.matrix(DNAmethylation_data),na.rm=TRUE))#37109 genes total in hgnc
  colnames(x)<-"Median_expr"
  DNAmethylation_data<-cbind(DNAmethylation_data,x)
  n<-ncol(DNAmethylation_data)
  n<-as.integer(n)
  fileNames1<-colnames(DNAmethylation_data)
  fileNames<-fileNames1[seq_len(n-1)]

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
    parallel::clusterExport(cl=cl,c("DNAmethylation_data"),envir=environment())#envir needed to be correct, export from global aswell as local
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
