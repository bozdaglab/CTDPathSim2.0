## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(CTDPathSim2.0)
library(printr)
data("CTDPathSim2")

## -----------------------------------------------------------------------------
knitr::kable(Methylation.Sample[1:5,1:3], digits = 2, caption = 'Probe-centric DNA methylation values of the tumor samples.')

## -----------------------------------------------------------------------------
knitr::kable(Reference.Methylation.Markers[1:5], caption = 'The marker loci (probes) of reference cell types.')

## -----------------------------------------------------------------------------
knitr::kable(Reference.Methylation.CellTypes[1:3,1:3], caption = 'R dataframe for the DNA methylation values for the reference cell types and the probes.')

## -----------------------------------------------------------------------------
knitr::kable(Annotation27K[1:3,], caption = 'An R dataframe object with probe, symbol and synonym columns of 27K DNA methylation data annotation.')

## -----------------------------------------------------------------------------
knitr::kable(Annotation450K[c(6,7,10),], caption = 'An R dataframe object with probe ID, UCSC_REFGENE_NAME and UCSC_REFGENE_GROUP columns 450K DNA methylation data annotation.')

## -----------------------------------------------------------------------------
knitr::kable(Expression.Sample[1:5,1:3], caption = 'An R dataframe for the gene-centric gene expression (RNASeq) values of the tumor samples.')

## -----------------------------------------------------------------------------
knitr::kable(Expression.CellLine[1:5,1:3], caption = 'An R dataframe for the gene-centric gene expression (RNASeq) values of the cell lines.')

## -----------------------------------------------------------------------------
knitr::kable(CancerGeneList[1:5,], caption = 'An R dataframe object with a column of frequently mutated cancer driver genes.')

## -----------------------------------------------------------------------------
knitr::kable(GeneCentric.DNAmethylation.Sample[1:5,1:3], caption = 'An R dataframe with gene-centric DNA methylation beta values of bulk tumor samples. Rows are the genes, columns are the tumor samples.')

## -----------------------------------------------------------------------------
knitr::kable(GeneCentric.DNAmethylation.CellLine[1:5,1:3], caption = 'An R dataframe with gene-centric DNA methylation beta values of cell lines. Rows are the genes, columns are the cell lines.')

## -----------------------------------------------------------------------------
knitr::kable(CNV.CellLine[1:5,1:3], caption = 'An R dataframe with gene-centric copy number values of cell lines. Rows are the genes, columns are the cell lines.')

## -----------------------------------------------------------------------------
?RunDeconvMethylFun

## ----results = "hide"---------------------------------------------------------
# Running the function
stage1_result_ct <-
  RunDeconvMethylFun(Methylation.Sample,Reference.Methylation.Markers)

## -----------------------------------------------------------------------------
# Outputs
knitr::kable(stage1_result_ct$methylation[1:5,], caption = 'Estimated DNA methylation profile of constituent cell types.')
knitr::kable(stage1_result_ct$proportions[1:5,], caption = 'Estimated proportions of constituent cell types.')

## -----------------------------------------------------------------------------
?PlotDeconvMethylFun

## -----------------------------------------------------------------------------
# Running the function
print(
PlotDeconvMethylFun(Methylation.Sample,Reference.Methylation.Markers,Reference.Methylation.CellTypes,stage1_result_ct,Reference.CellTypes.Names)
)
dev.off()

## -----------------------------------------------------------------------------
?ProbeToGeneFun
# Running the function
Deconv_methylation <- stage1_result_ct$methylation
Deconv_meth_gene <-
  ProbeToGeneFun(Annotation450K,Deconv_methylation,"450K")
# Output
knitr::kable(Deconv_meth_gene[1:5,], caption = 'An R dataframe object with deconvoluted DNA methylation values for each gene in cell types. The rows are the genes and the columns are the cell types. This function labels the columns as V1, V2, V3,....etc.')

## -----------------------------------------------------------------------------
?SampleMethylFun
# Running the function
CTDDirectory <- tempdir()
SampleMethylFun(Deconv_meth_gene,Deconv_proportions, CTDDirectory)
# Output
load("~/Sample_methylation/TCGA-OR-A5J1-01A.Rda")
knitr::kable(methylation[1:5,], caption = 'Deconvoluted DNA methylation profile of the sample TCGA-OR-A5J1-01A')

## -----------------------------------------------------------------------------
?RunDeconvExprFun
# Running the function
Deconv_expression<-RunDeconvExprFun(Expression.Sample,Deconv_proportions)
# Output
knitr::kable(Deconv_expression[1:5,], caption = 'An R dataframe object with deconvoluted gene expression values for each gene in cell types. The rows are the genes and the columns are the cell types. This function labels the columns as V1, V2, V3,....etc.')

## -----------------------------------------------------------------------------
?SampleExprFun
# Running the function
SampleExprFun(Deconv_expression,Deconv_proportions, CTDDirectory)
# Output
load("~/Sample_expression/TCGA-OR-A5J1-01A.Rda")
knitr::kable(expression[1:5,], caption = 'Deconvoluted gene expression profile of the sample TCGA-OR-A5J1-01A')

## -----------------------------------------------------------------------------
?GetSampleDEFun
# Running the function
CTDDirectory <- tempdir()
GetSampleDEFun(RnaSeq_data=Expression.Sample,parallel= TRUE,ncores=2, CTDDirectory)
# Output
load(paste0(CTDDirectory,"/Patient_DE_genes/DE_Gene/TCGA-OR-A5J6-01A_pat_de_genes.Rda"))
knitr::kable(pat_de_genes, caption = 'DE genes of the sample TCGA-OR-A5J1-01A')

## -----------------------------------------------------------------------------
?GetCellLineDEFun
# Running the function
GetCellLineDEFun(RnaSeq_data=Expression.CellLine,parallel= TRUE,ncores=2, CTDDirectory)
# Output
load(paste0(CTDDirectory,"/CellLine_DE_genes/DE_Gene/LK2_LUNG_cell_de_genes.Rda"))
knitr::kable(cell_de_genes[1:5,,drop=F], caption = 'DE genes of the LK2 lung cancer cell line')

## -----------------------------------------------------------------------------
?GetPathFun
# Running the function
data("CTDPathSim2")
Enriched.pathways.sample<-GetPathFun(CancerGeneList,pat_de_genes)
Enriched.pathways.cellLine<-GetPathFun(CancerGeneList,cell_de_genes)
# Output
knitr::kable(Enriched.pathways.sample[1:5,], caption = 'Sample(TCGA-OR-A5J1-01A)-specific enriched biological pathways')
knitr::kable(Enriched.pathways.cellLine[1:5,], caption = 'Cell line(LK2)-specific enriched biological pathways')


## -----------------------------------------------------------------------------
?GetSampleDMFun
# Running the function
GetSampleDMFun(GeneCentric.DNAmethylation.Sample,parallel= TRUE,ncores=2, CTDDirectory)
# Output
load(paste0(CTDDirectory,"/Patient_DM_genes/DM_Gene/TCGA-OR-A5J6-01A_pat_dm_genes.Rda"))
knitr::kable(pat_dm_genes, caption = 'DM genes of the sample TCGA-OR-A5J1-01A')

## -----------------------------------------------------------------------------
?GetCellLineDMFun
# Running the function
GetCellLineDMFun(GeneCentric.DNAmethylation.CellLine,parallel= TRUE,ncores=2, CTDDirectory)
# Output
load(paste0(CTDDirectory,"/CellLine_DM_genes/DM_Gene/LK2_LUNG_cell_dm_genes.Rda"))
knitr::kable(cell_dm_genes[1:5,], caption = 'DM genes of the LK2 cell line')

## -----------------------------------------------------------------------------
?GetCellLineDAFun
# Running the function
GetCellLineDAFun(CNV.CellLine,ncores=2, CTDDirectory)
# Output
load(paste0(CTDDirectory,"/CellLine_DA_genes/DA_Gene/LK2_LUNG_cell_dm_genes.Rda"))
knitr::kable(cell_da_genes[1:2,], caption = 'DA genes of the LK2 lung cancer cell line')

## -----------------------------------------------------------------------------
?FindSimFun
# Running the function to compute DNA methylation based Spearman similarity score of a patient-cell line pair
data("CTDPathSim2")
x<-FindSimFun(pat_dm_genes,pat_reactome,Deconv_meth_gene,cell_dm_genes,cell_reactome,ccle_methylation)
knitr::kable(x, caption = 'DNA methylation based Spearman similarity score of a patient-cell line pair')
# Running the function to compute gene expression-based Spearman similarity score of a patient-cell line pair
y<-FindSimFun(pat_de_genes,pat_reactome,Deconv_expression,cell_de_genes,cell_reactome,ccle_expr)
knitr::kable(y, caption = 'Gene expression-based Spearman similarity score of a patient-cell line pair')
# Running the function to compute copy number value based Spearman similarity score of a patient-cell line pair
z<-FindSimFun(pat_da_genes,pat_reactome,pat_cnv,cell_da_genes,cell_reactome,ccle_cnv)
knitr::kable(z, caption = 'Copy number aberration-based Spearman similarity score of a patient-cell line pair')

