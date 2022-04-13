#'Correlation plot to show the correlation between the reference cell types and estimated clusters DNA methylation deconvolution
#'
#'This function creates a correlation plot to show the correlation between the reference cell types used and estimated clusters after DNA methylation deconvolution. Users can use this plot to choose the appropriate labels for the cell types for the deconvoluted DNA methylation profiles.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 1 | Function 2 |

#' @param Methylation.Sample A R dataframe for the DNA methylation values for the samples and the probes.
#' @param Reference.Methylation.Probes A vector of the marker loci (probes) based on comparisons of each class of reference against all other samples using t-test.
#' @param Reference.Methylation.CellTypes A R dataframe for the DNA methylation values for the reference cell types and the probes.
#' @param stage1_result_ct The output from RunDeconvMethylFun()
#' @param Reference.CellTypes.Names A vector of the names of the reference cell types.
#'
#' @return Correlation Plot.
#' @export
#'
#' @examples
PlotDeconvMethylFun <-
  function(Methylation.Sample,
           Reference.Methylation.Probes,
           Reference.Methylation.CellTypes,
           stage1_result_ct,
           Reference.CellTypes.Names) {
    if ((is(Methylation.Sample, "SummarizedExperiment"))) {
      Methylation.Sample <- assay(Methylation.Sample)
    }
    if ((is(Reference.Methylation.CellTypes, "SummarizedExperiment"))) {
      Reference.Methylation.CellTypes <-
        assay(Reference.Methylation.CellTypes)
    }

    if ((is(Reference.CellTypes.Names, "SummarizedExperiment"))) {
      Reference.CellTypes.Names <- assay(Reference.CellTypes.Names)
    }


    chosen_Reference.Methylation.Probes <-
      intersect(row.names(Methylation.Sample),
                Reference.Methylation.Probes)
    # Compute correlation between estimated methylation profiles,
    # and reference methylation profiles
    cors_deconv_refs_ct = cor(Reference.Methylation.CellTypes[chosen_Reference.Methylation.Probes, ],
                              stage1_result_ct$methylation[chosen_Reference.Methylation.Probes, ])

    # Check what references had the highest correlation with each
    # of the estimated methylation profiles
    best_cors = rbind(apply(cors_deconv_refs_ct, 2, which.max),
                      apply(cors_deconv_refs_ct, 2, max))

    best_cor_labels = matrix("",
                             nrow = nrow(cors_deconv_refs_ct),
                             ncol = ncol(cors_deconv_refs_ct))

    for (i in seq_len(ncol(stage1_result_ct[["proportions"]]))) {
      best_cor_labels[best_cors[1, i], i] = as.character(round(best_cors[2, i], 2))
    }

    # Create a vector of colors representing the class of each reference
    ref_class_colors <-
      as.factor(as.character(Reference.CellTypes.Names$class))
    levels(ref_class_colors) <- RColorBrewer::brewer.pal(8, "Pastel1")
    ref_class_colors <- as.character(ref_class_colors)

    # Create a color gradient to be used in a heatmap of correlations
    color_gradient <- colorRampPalette(c("white", "steelblue"))
    # Plot correlation matrix
    par(mar=c(2,2,2,2))
    gplots::heatmap.2(
      cors_deconv_refs_ct,
      trace = "none",
      col = color_gradient(10),
      breaks = seq(0, 1, 0.1),
      margins = c(5, 5),
      RowSideColors = ref_class_colors,
      cellnote = best_cor_labels,
      notecol = "black"
    )

  }
