#' A data format to store genotype phenotype and sample QC results from SeqSQC.
#'
#' A SeqSQC object is a list of two objects. The first object \code{gdsfile} is the filepath of the GDS (discussed in section below) file which stores the genotype information from the original VCF file. The second object \code{QCresult} is a list of sample information and QC results, which include the dimension (# of samples and variants), sample annotation, and QC results for sample missing rate, sex check, inbreeding outlier check, IBD check, and population outlier check.
#'
#' @slot gdsfile A character string for the filepath of the GDS file. 
#' @slot QCresult A list with sample information and sample QC results. 
#' @name SeqSQC-class
#' @rdname SeqSQC-class
#' @aliases SeqSQC-class
#' @exportClass SeqSQC

## create class definitions
setClass("SeqSQC",
         ## contains = "gds.class",
         slots = c(
             gdsfile = "character",
             QCresult = "SimpleList")
         )

#' SeqSQC object Constructor
#' @name SeqSQC-class
#' @rdname SeqSQC-class
#' @param gdsfile A character string for the filepath of the GDS file.
#' @param QCresult A list with sample information and sample QC results.
#' @export SeqSQC
SeqSQC <- function(gdsfile, QCresult=List()){
    new("SeqSQC", gdsfile = gdsfile, QCresult = QCresult)
}

#' gdsfile getter and setter.
#' @rdname SeqSQC-class
#' @exportMethod gdsfile
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))

#' @rdname SeqSQC-class
#' @exportMethod "gdsfile<-"
setGeneric("gdsfile<-", function(x, ...) standardGeneric("gdsfile<-"))

#' QCresult getter and setter.
#' @rdname SeqSQC-class
#' @exportMethod QCresult
setGeneric("QCresult", function(x) standardGeneric("QCresult"))

#' @rdname SeqSQC-class
#' @exportMethod "QCresult<-"
setGeneric("QCresult<-", function(x, ...) standardGeneric("QCresult<-"))

#' @rdname SeqSQC-class
#' @aliases gdsfile,SeqSQC-method
#' @param x an SeqSQCClass object.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gdsfile(seqfile)
#' @return The filepath to the gds file. 
setMethod("gdsfile", "SeqSQC",function(x) x@gdsfile)

#' @rdname SeqSQC-class
#' @aliases "gdsfile<-",SeqSQC-method
setReplaceMethod("gdsfile", "SeqSQC", function(x, value) {
    x <- initialize(x, gdsfile = value)
    x
})  

#' @rdname SeqSQC-class
#' @name QCresult
#' @aliases QCresult,SeqSQC-method
#' @examples
#' QCresult(seqfile)
setMethod("QCresult", "SeqSQC", function(x) x@QCresult)

#' @rdname SeqSQC-class
#' @aliases "QCresult<-",SeqSQC-method
setReplaceMethod("QCresult", "SeqSQC", function(x, value) {
    x <- initialize(x, QCresult = value)
    x
})

## set show methods
setMethod("show", "SeqSQC",
          function(object){
              res <- object@QCresult
              cat("SeqSQC\n")
              cat("gds file:", object@gdsfile, "\n")
              cat("summary: 87 benchmark samples,", res$dimension[1]-87, "study samples,", res$dimension[2], "variants\n")
              cat("QC result:", paste(names(res), collapse=", "), "\n")
          }
          )

## test the validity of objects
setValidity("SeqSQC",
            function(object)
            {
                dat <- openfn.gds(object@gdsfile)
                
                if (!inherits(dat, "gds.class")){
                    return("object should inherit from 'gds.class'.")
                }
                var.names <- ls.gdsn(dat)
                sampleanno <- read.gdsn(index.gdsn(dat, "sample.annot"))
                pops <- unique(sampleanno$population)
                closefn.gds(dat)
                
                if (!all(c("sample.id", "sample.annot",
                           "snp.id", "snp.chromosome", "snp.position", "snp.allele", 
                           "genotype") %in% var.names)){
                    return("sample.id, sample.annot, snp.id, snp.chromosome, snp.position, snp.allele, and genotype are required variables.")
                }
                
                if(!all(c("sample", "population", "gender", "relation", "group") %in% names(sampleanno))){
                    return("sample annotation with sample name, population, gender, relation and group info are required")
                }
                
                if (!all(pops %in% c("AFR", "EUR", "EAS", "SAS", "ASN"))){
                    return("the population group must be AFR, EUR, EAS, SAS and ASN.")
                }
                TRUE
            }
            )


## load the benchmark data from ExperimentHub when installing the package
## write a fun <-  function() ExperimentHub::ExperimentHub()[["EH550"]]
## devtools::create("bigd")
## NAMESPACE: export("ExperimentHub")
 
## .onload <- function(libpath, pkgname){
##     message("libpath: ", libpath, " pkgname: ", pkgname)
##     ## hubdir <- file.path(libpath, pkgname, ".ExperimentHub")
##     ## supressMessages({
##     ## ExperimentHub::ExperimentHub(cache = hubdir)[["EH550"]]
##     ExperimentHub::ExperimentHub()[["EH550"]]
##     ## })
## }
## #' @ export

## myanalysis <- function(){
##     hubdir <- system.file(package="bigd", ".ExperimentHub")
##     mydata <- ExperimentHub::ExperimentHub(cache=hubdir)[["EH550"]]
##     message("benchmark data:", mydata$filename)
## }

