#'Enriched biological pathways
#'
#'This function identifies the enriched biological pathways for each sample and cell line utilizing the DE genes. It finds enriched biological pathways listed in the REACTOME database via a pathway enrichment tool Pathfinder. To prioritize cancer-related biological pathways, specifically, it uses the known frequently mutated genes that were considered cancer-driving genes from the DE gene list for this analysis. We obtained the frequently mutated cancer-driving genes from the COSMIC database. GetPath() considers the pathways that were significantly enriched with FDR-corrected hypergeometric p-values < 0.05 in this DE cancer gene list.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 3 | Function 3 |
#'
#' @param CancerGeneList A R dataframe with a column of frequently mutated cancer driver genes.
#' @param SampleCell_de_genes A R dataframe with a column containing DE genes of a sample or cell line.
#'
#' @return An R dataframe object with three columns labeled as 'ID', 'reactome_pathway', and 'p.adjust.pat', respectively.
#' @export
#'
#' @examples
GetPathFun<-function(CancerGeneList,SampleCell_de_genes){
  library(pathfindR)
  if((is(CancerGeneList, "SummarizedExperiment"))){
    CancerGeneList<-assay(CancerGeneList)
  }

  if((is(SampleCell_de_genes, "SummarizedExperiment"))){
    SampleCell_de_genes<-assay(SampleCell_de_genes)
  }

  common_gene_mut_de <-merge(SampleCell_de_genes,CancerGeneList,by.x = "gene",by.y = "V1")
  reactome<-enrichment(
    input_genes=common_gene_mut_de$gene,
    genes_by_term = pathfindR.data::reactome_genes,
    term_descriptions = pathfindR.data::reactome_descriptions,
    adj_method = "bonferroni",
    enrichment_threshold = 0.05,
    sig_genes_vec=common_gene_mut_de$gene,
    background_genes=unlist(pathfindR.data::reactome_genes))

  if (is.null(nrow(reactome))== TRUE){
    print("no patient reactome")
    pat_reactome<-data.frame("","","")
    colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
    return(pat_reactome)
  }else{
    d1<-data.frame(reactome$Term_Description)
    d2<-data.frame(reactome$adj_p)
    d3<-data.frame(reactome$ID)
    pat_reactome<-cbind(d3,d1,d2)
    colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
    return(pat_reactome)
  }
}


