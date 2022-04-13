#'Deconvolution of DNA methylation of bulk tumor samples
#'
#'To capture the accurate DNA methylation signal of tumor samples, this function utilizes a deconvolution-based algorithm to infer the samples deconvoluted DNA methylation profile with proportions of different cell types in the samples bulk tumor tissue. Deconvolution algorithm requires a reference DNA methylation profile of different cell types.Hence, we compiled a list of reference methylation profiles with eight different cell types, namely, B cell, natural killer, CD4T, CD8T, monocytes, generated from adult men with replicates.
#'
#'Created By: Banabithi Bose| Date Created: 4/10/2022 | Stage 1 | Function 1 |
#' @param Methylation.Sample A R dataframe for the probe-centric DNA methylation values of the tumor samples.
#' @param Reference.Methylation.Probes A vector of the marker loci (probes) based on comparisons of each class of reference against all other samples using t-test.
#'
#' @return A list of seven items: Estimated DNA methylation profile of constituent cell types, estimated proportions of constituent cell types, iteration number, explained variance, residual sum of squares, Akaike information criterion (AIC), and the residual sum of squares per iteration.
#' @export
#'
#' @examples
RunDeconvMethylFun <-
  function(Methylation.Sample,
           Reference.Methylation.Probes) {
    if ((is(Methylation.Sample, "SummarizedExperiment"))) {
      Methylation.Sample <- assay(Methylation.Sample)
    }

    chosen_Reference.Methylation.Probes <-
      intersect(row.names(Methylation.Sample),
                Reference.Methylation.Probes)

    stabilityResult <-
      estimate_stability(
        meth_bulk_samples = Methylation.Sample,
        informative_loci = chosen_Reference.Methylation.Probes,
        possible_num_ct = 3:8,
        subset_prop = 0.8,
        num_subsets = 10,
        reps_per_subset = 1,
        max_its = 800,
        rss_diff_stop = 1e-8
      )

    stage1_result_ct =
      run_edec_stage_1(
        meth_bulk_samples = Methylation.Sample,
        informative_loci = chosen_Reference.Methylation.Probes,
        num_cell_types = stabilityResult$most_stable_num_ct
      )

    stage1_result_ct$proportions <-
      round(stage1_result_ct$proportions, digits = 2)
    return(stage1_result_ct)
  }
