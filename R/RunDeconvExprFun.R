#'Deconvolution of gene expression of bulk tumor samples
#'
#'This function utilizes the estimated cell proportions from Step 1 as a fixed input and computes an average gene expression profiles of constituent cell types through a constrained least-squares using quadratic programming.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 2 | Function 1 |
#' @param Expression.Sample A R dataframe for the gene-centric gene expression (RNASeq) values of the tumor samples.
#' @param Deconv_proportions A R dataframe with estimated cell proportions in tumor samples from output of the RunDeconvMethylFun() function.
#'
#' @return A R dataframe with gene expression values in constituent cell types. The rows are the genes and the columns are the cell types. This function labels the columns as [,1], [,2],[,3],....etc.
#' @export
#'
#' @examples
RunDeconvExprFun<-function(Expression.Sample,Deconv_proportions){

  if((is(Expression.Sample, "SummarizedExperiment"))){
    Expression.Sample<-assay(Expression.Sample)
  }

  if((is(Deconv_proportions, "SummarizedExperiment"))){
    Deconv_proportions<-assay(Deconv_proportions)
  }

  x1<-data.table(colnames(Expression.Sample))
  x2<-data.table(row.names(Deconv_proportions))
  common_samples<-unique(merge(x1,x2,by="V1"))
  stage2_result_tumors = run_edec_stage_2(
    gene_exp_bulk_samples = Expression.Sample[,common_samples$V1],
    cell_type_props = Deconv_proportions[common_samples$V1,])
  Deconv_expression <-stage2_result_tumors$means
  return(Deconv_expression)
}
