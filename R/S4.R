## create class definitions
setClass("SeqSQCclass",
         ## contains = "gds.class",
         slots = c(
             gdsfile = "character",
             QCresult = "SimpleList")
         )

## Constructor
SeqSQCclass <- function(gdsfn, QCresult=List()){
    new("SeqSQCclass", gdsfile = gdsfn, QCresult = QCresult)
}

## set methods
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

