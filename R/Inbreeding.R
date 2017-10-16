#' Sample inbreeding check with SeqSQC object input file.
#'
#' Function to calculate population-specific inbreeding coefficients, and to predict inbreeding outliers that are five standard deviation beyond the mean. 
#' @param seqfile SeqSQC object, which includes the merged gds file for study cohort and benchmark.
#' @param remove.samples a vector of sample names for removal from inbreeding coefficient calculation. Could be problematic samples identified from previous QC steps, or user-defined samples.
#' @param LDprune whether to use LD-pruned snp set. The default is TRUE.
#' 
#' @param missing.rate to use the SNPs with "<= \code{missing.rate}" only; if NaN, no threshold. By default, we use \code{missing.rate = 0.1} to filter out variants with missing rate greater than 10\%.
#' @param ss.cutoff the minimum sample size (300 by default) to apply the MAF filter. This sample size is the sum of study samples and the benchmark samples of the same population as the study cohort.
#' @param maf to use the SNPs with ">= \code{maf}" if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise NaN is used by default for no MAF threshold.
#' @param hwe to use the SNPs with Hardy-Weinberg equilibrium p >= \code{hwe} if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise no hwe threshold. The default is 1e-6.
#' @param ... Arguments to be passed to other methods. 
#' @keywords inbreeding
#' @return a data frame with sample name, inbreeding coefficient, and an indicator of whether the inbreeding coefficient is five standard deviation beyond the mean. 
#' @details Using LD-pruned variants (by default), we calculate the inbreeding coefficients for each sample in the study cohort and for benchmark samples of the same population as the study cohort. Samples with inbreeding coefficients that are five standard deviations beyond the mean are considered problematic and are shown as "Yes" in the column of \code{outlier.5sd}. Benchmark samples in this column are set to be “NA”.
#' @export
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- Inbreeding(seqfile, remove.samples=NULL, LDprune=TRUE, missing.rate=0.1)
#' res.inb <- QCresult(seqfile)$Inbreeding
#' tail(res.inb)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

Inbreeding <- function(seqfile, remove.samples=NULL, LDprune=TRUE, missing.rate=0.1, ss.cutoff = 300, maf = 0.01, hwe = 1e-6, ...){

    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }
    message("calculating inbreeding coefficients ...")
    
    gfile <- SeqOpen(seqfile, readonly=TRUE)

    nds <- c("sample.id", "sample.annot", "snp.id") 
    allnds <- lapply(nds, function(x) read.gdsn(index.gdsn(gfile, x)))
    names(allnds) <- c("samples", "sampleanno", "snp.id")
    
    ## sample filters. (remove prespecified "remove.samples", and only keep samples within the study population)
    study.pop <- unique(allnds$sampleanno[allnds$sampleanno$group == "study", "population"])
    if(length(study.pop) > 1) stop("Study samples should be single population, please prepare input file accordingly.")
    
    if(study.pop == "ASN") study.pop <- c("EAS", "SAS", "ASN")
    
    if(!is.null(remove.samples)){
        flag <- allnds$sampleanno$population %in% study.pop & !allnds$samples %in% remove.samples
    }else{
        flag <- allnds$sampleanno$population %in% study.pop
    }
    sample.inb <- allnds$samples[flag]
    
    ## calculate population-specific inbreeding coefficient.

    ## add hwe filter for variants when sample size >= 300.
    if (length(sample.inb) >= ss.cutoff){
        snp.hwe <- snpgdsHWE(gfile, sample.id=sample.inb)
        hwe.idx <- snp.hwe > hwe & !is.na(snp.hwe)
    }else{
        hwe.idx <- rep(TRUE, length(allnds$snp.id))
    }
    
    ## use LDpruned SNPs if LDprune == TRUE.
    if(LDprune){
        ld <- read.gdsn(index.gdsn(gfile, "snp.annot/LDprune"))
    }else{
        ld <- rep(TRUE, length(allnds$snp.id))
    }
    
    ## SNP filters in together. (HWE+LDprune)
    snp.idx <- ld & hwe.idx
    
    ## use maf filter if sample size >= 300.
    if (length(sample.inb) >= ss.cutoff){ 
        rv <- snpgdsIndInb(gfile, snp.id=allnds$snp.id[snp.idx], sample.id=sample.inb, maf=maf, missing.rate=missing.rate, ...)
    }else{
        rv <- snpgdsIndInb(gfile, snp.id=allnds$snp.id[snp.idx], sample.id=sample.inb, maf=NaN, missing.rate=missing.rate, ...)
    }
    
    ## outliers
    m <- mean(rv$inbreeding)
    s <- sd(rv$inbreeding)
    coef <- 5
    outlier.5sd <- ifelse(rv$inbreeding < m-coef*s | rv$inbreeding > m+coef*s, "Yes", "No")
    
    res.inb <- data.frame(sample=sample.inb,  inbreeding=rv$inbreeding, outlier.5sd, stringsAsFactors=FALSE)

    ## set all outlier.* to be "NA" for 1000g benchmark samples. (We only use the 1000g samples to help identify inbreeding outliers from study cohort, and calculate inbreeding coef together with study samples within the sample population)
    skip.idx <- allnds$sampleanno[match(sample.inb, allnds$sampleanno$sample), 5] != "study"
    res.inb[skip.idx, "outlier.5sd"] <- NA
    closefn.gds(gfile)
    
    ## return the SeqSQC file with updated QC results.
    a <- QCresult(seqfile)
    a$Inbreeding <- res.inb
    outfile <- SeqSQC(gdsfile = gdsfile(seqfile), QCresult = a)
    return(outfile)
}
