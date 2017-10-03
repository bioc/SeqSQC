#' Sample missing rate check with SeqSQCclass input file.
#'
#' Function to calculate sample missing rate and to identify sample outlier with high missing rate (> 0.1).
#' @param seqfile SeqSQCclass input file, which includes the merged gds file for study cohort and benchmark.
#' @param remove.samples a vector of sample names for removal from missing rate check. Could be problematic samples identified from other QC steps, or user-defined samples.
#' @keywords MissingRate
#' @return a data frame with sample name, sample missing rate, and an indicator of whether the sample has a missing rate greater than 0.1.
#' @details The value of the outlier column is set to NA for benchmark samples.
#' @export
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQCclass(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- MissingRate(seqfile, remove.samples=NULL)
#' res.mr <- QCresult(seqfile)$MissingRate
#' tail(res.mr)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

MissingRate <- function(seqfile, remove.samples=NULL){

    ## check
    if (!inherits(seqfile, "SeqSQCclass")){
        return("object should inherit from 'SeqSQCclass'.")
    }

    message("calculating sample missing rates...")
    
    gfile <- SeqOpen(seqfile, readonly=TRUE)
    samples <- read.gdsn(index.gdsn(gfile, "sample.id"))
    sampleanno <- read.gdsn(index.gdsn(gfile, "sample.annot"))

    ## sample filters. (remove prespecified "remove.samples", and only keep samples within the study population)
    study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
    if(length(study.pop) > 1) stop("Study samples should be single population, please prepare input file accordingly.")
    if(study.pop == "ASN") study.pop <- c("EAS", "SAS", "ASN")
    
    if(!is.null(remove.samples)){
        flag <- sampleanno$population %in% study.pop & !samples %in% remove.samples
    }else{
        flag <- sampleanno$population %in% study.pop
    }
    sample.mr <- samples[flag]
    
    gt <- readex.gdsn(index.gdsn(gfile, "genotype"), list(NULL, samples %in% sample.mr))
    mr <- apply(gt, 2, function(x) sum(! x %in% c(0, 1, 2))/length(x))
    outlier <- ifelse(mr > 0.1, "Yes", "No")
    res.mr <- data.frame(sample = sample.mr, missingRate = mr, outlier = outlier, stringsAsFactors = FALSE)

    skip.idx <- sampleanno[match(sample.mr, sampleanno$sample), 5] != "study"
    res.mr[skip.idx, "outlier"] <- NA

    closefn.gds(gfile)

    ## return the SeqSQCclass file with updated QC results.
    a <- QCresult(seqfile)
    a$MissingRate <- res.mr
    ## QCresult(seqfile)$MissingRate <-  res.mr
    outfile <- SeqSQCclass(gdsfile = gdsfile(seqfile), QCresult = a)
    return(outfile)
}

