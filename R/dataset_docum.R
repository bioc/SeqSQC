#' example gds file used in vignette.
#'
#' This gds file contains genotype and phenotype for 92 whole-genome sequenced samples captured by CCDS region.
#' This is a merged dataset of the 87 benchmark samples and the 5 study samples (all are assembled from the 1000 Genomes Project).
#' The meta info for these 92 samples includes sample name, pupulation, age, relation note and group info (benchmark or study). 
#'
#' @name example.gds
#' @docType data
#' @author Qian Liu \email{qliu7@@buffalo.edu}
#' @keywords datasets
NULL

#' Example vcf file used in vignette.
#'
#' This vcf file contains only a subset (1000 lines of variants) of the original vcf file for the 5 study samples (examples assembled from the 1000 Genomes Project).
#' This is to be used as a runnable example in the function of \code{LoadVfile} and \code{sampleQC} in the vignette. 
#'
#' @name example_sub.vcf
#' @docType data
#' @author Qian Liu \email{qliu7@@buffalo.edu}
#' @keywords datasets
NULL

#' Example SeqSQC file used in vignette.
#'
#' The SeqSQC object is a list of two objects.
#' The first object \code{gdsfile} is the filepath of the "example.gds" file which stores the genotype and meta info of the example data merged with the benchmark data.
#' The second object \code{QCresult} contains the data dimensions (# of samples and variants), sample annotation, and QC results for sample missing rate, sex check, inbreeding outlier check, IBD check, and population outlier check.
#' @name example.seqfile.Rdata
#' @docType data
#' @author Qian Liu \email{qliu7@@buffalo.edu}
#' @keywords datasets
NULL

#' Sample annotation file for the example data used in vignette.
#'
#' This sample annotation file is a required input from the user when using SeqSQC.
#' It includes the sample info with sample name stored in the column of \code{sample}, the population info stored in the column of \code{population}, and the gender info stored in the column of \code{gender}.
#' The \code{population} column must be a in the format of "AFR/EUR/ASN/EAS/SAS". The \code{gender} column must be in the format of "female/male". 
#' @name sampleAnnotation.txt
#' @docType data
#' @author Qian Liu \email{qliu7@@buffalo.edu}
#' @keywords datasets
NULL
