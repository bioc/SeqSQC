#' A data format to store genotype phenotype and sample QC results from SeqSQC.
#'
#' A SeqSQCclass object is a list of two objects. The first object \code{gdsfile} is the filepath of the GDS (discussed in section below) file which stores the genotype information from the original VCF file. The second object \code{QCresult} is a list of sample information and QC results, which include the dimension (# of samples and variants), sample annotation, and QC results for sample missing rate, sex check, inbreeding outlier check, IBD check, and population outlier check.
#'
#' @slot gdsfile A character string for the filepath of the GDS file. 
#' @slot QCresult A list with sample information and sample QC results. 
#' @name SeqSQCclass-class
#' @rdname SeqSQCclass-class
#' @aliases SeqSQCclass-class
#' @exportClass SeqSQCclass

## create class definitions
setClass("SeqSQCclass",
         ## contains = "gds.class",
         slots = c(
             gdsfile = "character",
             QCresult = "SimpleList")
         )

#' SeqSQCclass object Constructor
#' @name SeqSQCclass-class
#' @rdname SeqSQCclass-class
#' @param gdsfile A character string for the filepath of the GDS file.
#' @param QCresult A list with sample information and sample QC results.
#' @export SeqSQCclass
SeqSQCclass <- function(gdsfile, QCresult=List()){
    new("SeqSQCclass", gdsfile = gdsfile, QCresult = QCresult)
}

#' Method gdsfile.
#' @name SeqSQCclass-class
#' @rdname SeqSQCclass-class
#' @exportMethod gdsfile
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))

#' Method QCresult.
#' @name SeqSQCclass-class
#' @rdname SeqSQCclass-class
#' @exportMethod QCresult
setGeneric("QCresult", function(x) standardGeneric("QCresult"))

#' @rdname SeqSQCclass-class
#' @name gdsfile
#' @aliases gdsfile,SeqSQCclass-method
#' @param x an SeqSQCClass object.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gdsfile(seqfile)
#' @return The filepath to the gds file. 

setMethod("gdsfile", "SeqSQCclass",function(x) x@gdsfile)

#' @rdname SeqSQCclass-class
#' @name QCresult
#' @aliases QCresult,SeqSQCclass-method
#' @examples
#' QCresult(seqfile)
setMethod("QCresult", "SeqSQCclass", function(x) x@QCresult)

## #' @rdname gdsfile
## #' @export
## setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))

## #' @rdname QCresult
## #' @export
## setGeneric("QCresult", function(x) standardGeneric("QCresult"))

## #' Accessors for the 'gdsfile' slot of a seqSQCclass object.
## #' the gdsfile slot holds the filepath of the gds file, which contains the genotypes and meta info for samples and variants.
## #' @docType methods
## #' @param x a SeqSQCclass object
## #' @name gdsfile
## #' @rdname gdsfile
## #' @aliases gdsfile,SeqSQCclass-method
## #' @exportMethod gdsfile
## setMethod("gdsfile", "SeqSQCclass",
##           function(x){
##               gp <- x@gdsfile
##               return(gp)
##           }
##           )

## #' Accessors for the 'QCresult' slot of a SeqSQCclass object.
## #' the QCresult slot holds the sample QC results, including a summary of variants and samples, and QC results from missing rate check, sex check, inbreeding outlier check, IBD check and population outlier check.
## #' @docType methods
## #' @param x a SeqSQCclass object
## #' @name QCresult
## #' @rdname QCresult
## #' @aliases QCresult,SeqSQCclass-method
## #' @exportMethod QCresult
## setMethod("QCresult", "SeqSQCclass",
##           function(x){
##               res <- x@QCresult
##               return(res)
##           }
##           )

## set show methods
setMethod("show", "SeqSQCclass",
          function(object){
              res <- object@QCresult
              cat("SeqSQCclass\n")
              cat("gds file:", object@gdsfile, "\n")
              cat("summary: 87 benchmark samples,", res$dimension[1]-87, "study samples,", res$dimension[2], "variants\n")
              cat("QC result:", paste(names(res), collapse=", "), "\n")
          }
          )

## test the validity of objects
setValidity("SeqSQCclass",
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

## todo: update QC steps, finished and return a results, not gds file. 
## todo: update sampleQC.R, gfile = "" (character); close the gds file when QC. Open inside QC.
## todo: 

## set methods
## setGeneric("plotMissingRate", function(x) standardGeneric("plotMissingRate"))
## setMethod("plotMissingRate", "sampleQCGDSClass", function() )

## Defining a coercion method
## setAs("SeqSQCclass", "list",
##         function(from)
##         {
##             nodes <- ls.gdsn(index.gdsn(from, "results"))
##             res <- list()
##             for (i in 1:length(nodes)){
##                 res[[i]] <- read.gdsn(index.gdsn(from, paste0("results/", nodes[i])))
##             }
##             names(res) <- nodes
##             res
##         }
##       )
## results <- as(dat, "list") # testing the gds object "dat" (which is opened in R console)

