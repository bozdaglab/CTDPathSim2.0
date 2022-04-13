#'Cell line-specific differentially aberrated (DA) genes
#'
#'This function selects the highly amplified and deleted genes for each cell line compared with all other cell lines from various cancer types to capture the cancer-specific variation in copy number profiles in computing the similarity score between a patient sample and cell line. It selects the cell line-specific DA genes based on Tukey’s mean-difference curve based on gene-centric copy number data of cell lines. To determine whether a gene is DA for a specific cell line, it calculates the mean of that gene's copy number values across all cell lines and a difference between the copy number value of that gene in the cell line under consideration and the mean of the copy number values of that gene in the other cell lines. For each cell line, the highly amplified genes are chosen based on mean ≥ 0.5 and difference ≤ 1, and highly deleted genes are chosen based on mean ≥ 0.5 and difference ≤ -1.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 4 | Function 1 |
#' @param CNV.CellLine A R dataframe with gene-centric copy number values of cell lines. Rows are the genes, columns are the cell lines.
#' @param ncores An integer value for the number of cores.
#' @param CTDDirectory A character string for the file path of the directory for the output to be stored.
#'
#' @return R dataframe objects each with a single column containing DA genes for each cell line. The column is labeled as gene.
#' @export
#'
#' @examples
GetCellLineDAFun <- function(CNV.CellLine,
                             ncores = 1,
                             CTDDirectory = "~") {
  if ((is(CNV.CellLine, "SummarizedExperiment"))) {
    CNV.CellLine <- assay(CNV.CellLine)
  }
  if (!file.exists(file.path(CTDDirectory, "CellLine_DA_genes"))) {
    dir.create(file.path(CTDDirectory, "CellLine_DA_genes"))
  }
  if (!file.exists(file.path(CTDDirectory, "CellLine_DA_genes", "UP_Gene"))) {
    dir.create(file.path(CTDDirectory, "CellLine_DA_genes", "UP_Gene"))
  }
  if (!file.exists(file.path(CTDDirectory, "CellLine_DA_genes", "DOWN_Gene"))) {
    dir.create(file.path(CTDDirectory, "CellLine_DA_genes", "DOWN_Gene"))
  }
  if (!file.exists(file.path(CTDDirectory, "CellLine_DA_genes", "DA_Gene"))) {
    dir.create(file.path(CTDDirectory, "CellLine_DA_genes", "DA_Gene"))
  }

  processFiles <- function(i) {
    d <-
      data.frame(CNV.CellLine[, i] - rowMeans(CNV.CellLine[,-i]))
    m <- data.frame(rowMeans(CNV.CellLine[,]))
    df <- cbind(d, m)
    df <- na.omit(df)
    colnames(df) <- c("d", "m")
    df$z <- "ok"
    df[df$d >= 1 & df$m >= .5, 3] <- "upda"
    df[df$d <= (-1) & df$m >= .5, 3] <- "downda"
    SampleName <- colnames(CNV.CellLine[, i, drop = FALSE])
    ampgene <- data.table(row.names(df[df$z == "upda",]))
    delgene <- data.table(row.names(df[df$z == "downda",]))
    cell_da_genes <- data.table(rbind(ampgene, delgene))

    save(ampgene,
         file = file.path(
           CTDDirectory,
           "CellLine_DA_genes",
           "UP_Gene",
           paste0(SampleName, "_upgene.Rda")
         ))
    save(delgene,
         file = file.path(
           CTDDirectory,
           "CellLine_DA_genes",
           "DOWN_Gene",
           paste0(SampleName, "_downgene.Rda")
         ))
    save(cell_da_genes,
         file = file.path(
           CTDDirectory,
           "CellLine_DA_genes",
           "DA_Gene",
           paste0(SampleName, "_cell_dm_genes.Rda")
         ))
  }
  n <- ncol(CNV.CellLine)
  n <- as.integer(n)
  i <- as.integer(1:n)
  if (ncores > 1) {
    print("cluster making started")
    cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
    parallel::clusterExport(cl = cl, NULL, envir = environment())
    parLapplyLB(cl, i, processFiles)
    stopCluster(cl)

  } else {
    lapply(i, processFiles)
  }
  return()
}
