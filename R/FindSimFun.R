#'Sample-cell line pathway activity-based similarity score
#'
#'This function computes Spearman rank correlation between each sample-cell line pair using the sample-specific deconvoluted expression, sample-specific deconvoluted DNA methylation, and gene-centric copy number values of samples. All DM genes that occur in an enriched pathway of a sample or cell line are used to compute the Spearman rank correlation between each sample-cell line pair. Similarly, using DE and DA genes that are present in the gene list of the union of enriched pathways, it computes the Spearman rank correlation between each sample-cell line pair. Using DM and DE genes, it uses the deconvoluted profiles of samples to compute Spearman rank correlation whereas for DA genes, it uses the gene-centric copy number profile of samples to compute Spearman rank correlation.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 5 | Function 1 |
#' @param Patient.D.Genes A vector of DE or DM or DA genes of a tumor sample
#' @param Patient.Reactome A R dataframe for patient-specific enriched biological pathways with three columns labeled as ID, reactome_pathway, and p.adjust.pat, respectively.
#' @param Patient.Expr.Meth.CNV A R dataframe with sample-specific deconvoluted gene expression or deconvoluted DNA methylation or copy number profiles.
#' @param Cell.D.Genes A vector of DE or DM or DA genes of a cell line
#' @param Cell.Reactome A R dataframe for cell line-specific enriched biological pathways with three columns labeled as ID, reactome_pathway, and p.adjust.pat, respectively.
#' @param Cell.Expr.Meth.CNV A R dataframe with cell line-specific gene expression or DNA methylation or copy number profiles.
#'
#' @return A numeric value (Spearman similarity score)
#' @export
#'
#' @examples
FindSimFun<- function(Patient.D.Genes,Patient.Reactome,Patient.Expr.Meth.CNV,Cell.D.Genes,Cell.Reactome,Cell.Expr.Meth.CNV){

  if((is(Patient.D.Genes, "SummarizedExperiment"))){
    Patient.D.Genes<-assay(Patient.D.Genes)
  }

  if((is(Patient.Reactome, "SummarizedExperiment"))){
    Patient.Reactome<-assay(Patient.Reactome)
  }

  if((is(Patient.Expr.Meth.CNV, "SummarizedExperiment"))){
    Patient.Expr.Meth.CNV<-assay(Patient.Expr.Meth.CNV)
  }

  if((is(Cell.D.Genes, "SummarizedExperiment"))){
    Cell.D.Genes<-assay(Cell.D.Genes)
  }

  if((is(Cell.Reactome, "SummarizedExperiment"))){
    Cell.Reactome<-assay(Cell.Reactome)
  }

  if((is(Cell.Expr.Meth.CNV, "SummarizedExperiment"))){
    Cell.Expr.Meth.CNV<-assay(Cell.Expr.Meth.CNV)
  }

  SampleName<-NULL
  cellLine_id<-NULL
  tryCatch(
    {
      row_pat<-row.names(Patient.Expr.Meth.CNV)
      row_cell<-row.names(Cell.Expr.Meth.CNV)
      row_pat_cell<-intersect(row_pat,row_cell)##common genes
      union_reactome<-union(Patient.Reactome$ID,Cell.Reactome$ID)## union of pathways bet patient-cell line
      patients_DE<-intersect(Patient.D.Genes$gene,row_pat_cell)
      cellLines_DE<-intersect(Cell.D.Genes$gene,row_pat_cell)
      union_de_genes1<-union(Patient.D.Genes$gene,Cell.D.Genes$gene)# union DE genes bet patient cell line
      union_de_genes<-intersect(union_de_genes1,row_pat_cell)

      Lx<-Patient.Expr.Meth.CNV[union_de_genes,,drop=F]## patient's Patient.Expr.Meth.CNV with common DE genes
      r1<-data.table(row.names(Lx))
      colnames(r1)<-"g"
      Lx<-cbind(r1,Lx)
      Ly<-as.data.frame(Cell.Expr.Meth.CNV[union_de_genes,,drop=FALSE])## cell line's Patient.Expr.Meth.CNV with union DE genes
      r2<-data.table(row.names(Ly))
      colnames(r2)<-"g"
      Ly<-cbind(r2,Ly)
      ### For reactome
      reactome_pathNames<-union_reactome
      if(length(reactome_pathNames)==0){
        reactome_pathNames<-"NA"
        reactome_descriptions<-"NA"
      }else{
        reactome_genes_by_term = data.frame(cbind(pathfindR.data::reactome_genes))
        reactome_term_descriptions = data.frame(pathfindR.data::reactome_descriptions)
        reactome_descriptions <- as.character(reactome_term_descriptions[reactome_pathNames,])
        reactome_pathway_genes<-data.table(unique(unlist(reactome_genes_by_term[reactome_pathNames,])))## taking union of genes in all intersected pathways
        colnames(reactome_pathway_genes)<-"g"
        Common_DE_genes_pathway <-unique(intersect(reactome_pathway_genes$g,union_de_genes))
      }
      if (nrow(reactome_pathway_genes)==0){
        reactomeLx<-matrix(1,nrow = 0,ncol=1)
      }else{

        colnames(reactome_pathway_genes)<-"g"
        reactomeLx <- merge(reactome_pathway_genes,Lx,by="g")
        reactomeLy <- merge(reactome_pathway_genes,Ly,by="g")
      }
      if (nrow(reactomeLx)==0||nrow(reactomeLx)==1){

        print("corr not possible")

        reactome_spearScore<-"NA"

        Similarity_score<-cbind(paste0(str_sub(SampleName,1,-5)),cellLine_id,reactome_pathNames,reactome_descriptions,reactome_spearScore,length(patients_DE),length(cellLines_DE),
                                        length(unique(union_de_genes)),length(unique(reactome_pathway_genes$g)),length(unique(Common_DE_genes_pathway)))
        Similarity_score<-data.table(Similarity_score)
        colnames(Similarity_score)<-c("patient","cellLine","reactome_pathway","description","Spearman_Similarity_Score","pat_de","cell_de","union_de","reactome_genes","de_in_reactome")
        Similarity_score<-unique(Similarity_score[,c(1,2,5:10)])

      }
      ## For CNV
      if (ncol(Patient.Expr.Meth.CNV)==1){
        reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
        x1<-data.frame(reactome_spearCorr1[1])
        colnames(x1)<-"V1"
      }
      if (ncol(Patient.Expr.Meth.CNV)==3){
        reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
        reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
        reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
        x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1]))
        colnames(x1)<-"V1"
      }
      if (ncol(Patient.Expr.Meth.CNV)==4){
        reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
        reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
        reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
        reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
        x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1],
                             reactome_spearCorr4[1]))
        colnames(x1)<-"V1"
      }
      if (ncol(Patient.Expr.Meth.CNV)==5){
        reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
        reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
        reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
        reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
        reactome_spearCorr5<-cor(reactomeLy[,2],reactomeLx[,6],use ="complete",method=c("spearman"))
        x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1]),reactome_spearCorr4[1],reactome_spearCorr5[1])
        colnames(x1)<-"V1"
      }

      if(abs(mean(x1$V1,na.rm=TRUE))=="NaN"){
        reactome_spearScore<-"NA"
      }else{
        reactome_spearScore<-abs(mean(x1$V1,na.rm=TRUE))
      }

      Similarity_score<-cbind(reactome_spearScore,length(patients_DE),length(cellLines_DE),
                                      length(unique(union_de_genes)),length(unique(reactome_pathway_genes$g)),length(unique(Common_DE_genes_pathway)))
      Similarity_score<-data.table(Similarity_score)
      colnames(Similarity_score)<-c("Spearman_Similarity_Score","pat_de","cell_de","union_de","reactome_genes","de_in_reactome")
    },
    error = function(error_condition) {
    return(Similarity_score[,1,drop=F])
    },finally={
      return(Similarity_score[,1,drop=F])
    })
}
