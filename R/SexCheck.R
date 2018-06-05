#' Sample gender check with SeqSQC object input file.
#'
#' Function to calculate the X chromosome inbreeding coefficient and
#' to predict sample gender.
#' @param seqfile SeqSQC object, which includes the merged gds file
#'     for study cohort and benchmark.
#' @param remove.samples a vector of sample names for removal from sex
#'     check. Could be problematic samples identified from previous QC
#'     steps, or user-defined samples.
#' @param missing.rate to use the SNPs with "<= \code{missing.rate}"
#'     only; if NaN, no threshold. By default, we use
#'     \code{missing.rate = 0.1} to filter out variants with missing
#'     rate greater than 10\%.
#' @param ss.cutoff the minimum sample size (300 by default) to apply
#'     the MAF filter. This sample size is the sum of study samples
#'     and the benchmark samples of the same population as the study
#'     cohort.
#' @param maf to use the SNPs with ">= \code{maf}" if sample size
#'     defined in above argument is greater than \code{ss.cutoff};
#'     otherwise NaN is used by default for no MAF threshold.
#' @param ... Arguments to be passed to other methods.
#' @keywords SexCheck
#' @return a data frame with sample name, reported gender, x
#'     chromosome inbreeding coefficient, and predicted gender.
#' @details Samples are predicted to be female or male if the
#'     inbreeding coefficient is below 0.2, or greater than 0.8,
#'     respectively. The samples with discordant reported gender and
#'     predicted gender are considered as problematic. When the
#'     inbreeding coefficient is within the range of [0.2, 0.8], “0”
#'     is shown in the column of \code{pred.sex} to indicate ambiguous
#'     gender, which is not considered as problematic.
#' @export
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- SexCheck(seqfile, remove.samples=NULL, missing.rate=0.1)
#' res.sexc <- QCresult(seqfile)$SexCheck
#' tail(res.sexc)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

SexCheck <- function(seqfile, remove.samples=NULL, missing.rate = 0.1,
                     ss.cutoff = 300, maf = 0.01, ...){

    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }

    message("calculating sex inbreeding ...")
    
    gfile <- SeqOpen(seqfile, readonly=TRUE)
    
    nds <- c("sample.id", "sample.annot", "snp.chromosome", "snp.position", "snp.id") 
    allnds <- lapply(nds, function(x) read.gdsn(index.gdsn(gfile, x)))
    names(allnds) <- c("samples", "sampleanno", "snp.chr", "snp.pos", "snps")

     ## sample filters. (remove prespecified "remove.samples", and only keep samples within the study population)
    study.pop <- unique(allnds$sampleanno[allnds$sampleanno$group == "study", "population"])
    if (length(study.pop) > 1)
        stop("Study samples should be single population, please prepare input file accordingly.")

    if(study.pop == "ASN") study.pop <- c("EAS", "SAS", "ASN")
    
    if(!is.null(remove.samples)){
        flag <- allnds$sampleanno$population %in% study.pop & !allnds$samples %in% remove.samples
    }else{
        flag <- allnds$sampleanno$population %in% study.pop
    }
    sample.sex <- allnds$samples[flag]
        
    ## SNP filters. (only keep X variants, and remove X variants in pseudo-autosomal region.)
    flag <- allnds$snp.chr == "X"
    snp.x <- allnds$snps[flag]
    snp.chr.x <- allnds$snp.chr[flag]
    snp.pos.x <- allnds$snp.pos[flag]

    ## remove pseudo-autosomal region in X chromosome.
    grx <- GRanges(snp.chr.x, IRanges(snp.pos.x, snp.pos.x))
    par <- GRanges(c("X", "X"), IRanges(start=c(60001, 154931044), end=c(2699520, 155260560)))
    ov <- findOverlaps(grx, par)
    par.rm <- unique(queryHits(ov))
    snp.x.final <- snp.x[-par.rm]

    ## calculate population-specific sex inbreeding coefficient.

    ## use maf filter for variants when sample size >= 300.
    if (length(sample.sex) >= ss.cutoff){
        sexinb <- snpgdsIndInb(gfile, autosome.only=FALSE,
                               sample.id=sample.sex,
                               snp.id=snp.x.final, maf=maf,
                               missing.rate=missing.rate, ...)
    }else{
        sexinb <- snpgdsIndInb(gfile, autosome.only=FALSE,
                               sample.id=sample.sex,
                               snp.id=snp.x.final, maf=NaN,
                               missing.rate=missing.rate, ...)
    }
    closefn.gds(gfile)
    
    ## extract gender info
    sex <- allnds$sampleanno[match(sample.sex, allnds$sampleanno[,1]), 3]
    res.sexcheck <- data.frame(sample=sample.sex, sex=sex,
                               sexinb=sexinb$inbreeding,
                               stringsAsFactors=FALSE)
    
    ## gender prediction.
    pred <- ifelse(res.sexcheck$sexinb > 0.8, "male",
            ifelse(res.sexcheck$sexinb < 0.2, "female", 0))
    res.sexcheck$pred.sex <- pred
    
    ## return the SeqSQC file with updated QC results.
    a <- QCresult(seqfile)
    a$SexCheck <- res.sexcheck
    outfile <- SeqSQC(gdsfile = gdsfile(seqfile), QCresult = a)
    return(outfile)
}
