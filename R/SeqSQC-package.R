#' SeqSQC
#'
#' Sample Quality Check for NGS Data.
#'
#' @name SeqSQC-package
#' @title Sample Quality Check for NGS Data using SeqSQC package
#' @import gdsfmt
#' @import SNPRelate
#' @import ExperimentHub
#' @author Qian Liu
#' @seealso
#' \code{\link{LoadVfile}} \cr
#' for data preparation; \cr
#' \code{\link{MissingRate}} \cr
#' \code{\link{PCACheck}} \cr
#' \code{\link{Inbreeding}} \cr
#' \code{\link{IBDCheck}} \cr
#' \code{\link{PCACheck}} \cr
#' for individual sample QC checks; \cr
#' \code{\link{problemList}} \cr
#' for the summary of problematic samples with reason and sample list to be removed; \cr
#' \code{\link{IBDRemove}} \cr
#' for the problematic sample pairs detected with cryptic relationship; \cr
#' \code{\link{RenderReport}} \cr
#' to generate the sample QC report; \cr
#' \code{\link{plotQC}} \cr
#' to generate the ggplot or interactive plots in html format for each individual QC check; \cr
#' \code{\link{sampleQC}} \cr
#' for wrapper of data preparation, all sample QC checks, QC result summary, and sample QC report.  

NULL
