#'Converting DNA methylation data from probes to genes
#'
#'This function converts the DNA methylation data from probe-centric to gene-centric. This function works for both the 450K probe format and 27K probe format.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 1 | Function 3 |
#' @param Methylation.Probe.Annotation A R dataframe object with 'probe', 'symbol' and 'synonym' columns for 27K or 'ID', 'UCSC_REFGENE_NAME', 'UCSC_REFGENE_GROUP' columns for 450K DNA methylation data according to the convention of the annotation file in "Genome Browser" (https://genome.ucsc.edu/.
#' @param Deconv_methylation A R dataframe with DNA methylation values for each probe and cell type.
#' @param Annotation.skim A character string indicating if the user is using either 450K or 27K probes.
#'
#' @return An R dataframe object with deconvoluted DNA methylation values for each gene in cell types. The rows are the genes and the columns are the cell types. This function labels the columns as V1, V2, V3,....etc.
#' @export
#'
#' @examples
ProbeToGeneFun <-
  function(Methylation.Probe.Annotation,
           Deconv_methylation,
           Annotation.skim) {
    library(data.table)
    library(tidyr)
    library(sqldf)
    if ((is(Methylation.Probe.Annotation, "SummarizedExpirament"))) {
      Methylation.Probe.Annotation <- assay(Methylation.Probe.Annotation)
    }

    if ((is(Deconv_methylation, "SummarizedExpirament"))) {
      Deconv_methylation <- assay(Deconv_methylation)
    }

    if (Annotation.skim != "27K" & Annotation.skim != "450K") {
      message("Incorrect gene format input. Input either 27K or 450K.")
    }


    if (Annotation.skim == "27K") {
      A_27 <- Methylation.Probe.Annotation
      OV_Deconv_methylation <- NULL
      ###############Process Beta Values for 27K###############
      probes <- data.table(row.names(Deconv_methylation))
      colnames(probes) <- "ID"
    } else if (Annotation.skim == "450K") {
      A_450 <- Methylation.Probe.Annotation

      An1 <-
        sqldf::sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS1500%' ")#59310
      An2 <-
        sqldf::sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS200%' ")#46466
      An3 <-
        sqldf::sqldf("Select * from A_450 where UCSC_RefGene_Group like '%5UTR%' ")#38876
      Probe_Anno_450K <- rbind(An1, An2, An3)#144652
    }
    probes <- data.table(row.names(Deconv_methylation))
    colnames(probes) <- "ID"
    Deconv_methylation <- data.table(cbind(probes, Deconv_methylation))

    if (Annotation.skim == "27K") {
      B <- merge(A_27, Deconv_methylation, by = "ID")
      B$gene <- paste(B$Symbol, ";", B$Synonym)
      B$D <-
        vapply(strsplit(as.character(B$gene), ";", fixed = TRUE), function(Methylation.Probe.Annotation)
          paste(unique(Methylation.Probe.Annotation), collapse = ";"), FUN.VALUE =
            character(1))
    } else if (Annotation.skim == "450K") {
      B <- merge(Probe_Anno_450K, Deconv_methylation, by = "ID")#125139
      B$D <-
        vapply(strsplit(as.character(B$UCSC_REFGENE_NAME), ";", fixed = TRUE), function(Methylation.Probe.Annotation)
          paste(unique(Methylation.Probe.Annotation), collapse = ";"), FUN.VALUE =
            character(1))
    }

    D_Probe.Annotation <- separate_rows(B, D, convert = TRUE)
    D_Probe.Annotation <- data.frame(D_Probe.Annotation)
    n1 <- ncol(D_Probe.Annotation)
    if (Annotation.skim == "27K") {
      n2 <- n1 - 2
    } else if (Annotation.skim == "450K") {
      n2 <- n1 - 1
    }
    D_Probe.Annotation_1 <- D_Probe.Annotation[, c(n1, 4:n2)]
    Deconv_meth_gene <-
      aggregate(
        D_Probe.Annotation_1[,-1],
        by = list(D_Probe.Annotation_1$D),
        mean,
        na.rm = TRUE
      )

    row.names(Deconv_meth_gene) <- Deconv_meth_gene$Group.1
    Deconv_meth_gene <- Deconv_meth_gene[, -1]

    return(Deconv_meth_gene)



  }
