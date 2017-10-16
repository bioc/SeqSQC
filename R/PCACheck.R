#' Population outlier check with SeqSQC object input file.
#'
#' Function to perform principle component analysis for all samples and to infer sample ancestry.
#' @param seqfile SeqSQC object, which includes the merged gds file for study cohort and benchmark.
#' @param remove.samples a vector of sample names for removal from PCA calculation. Could be problematic samples identified from previous QC steps, or user-defined samples.
#' @param LDprune whether to use LD-pruned snp set, the default is TRUE.
#' @param missing.rate to use the SNPs with "<= \code{missing.rate}" only; if NaN, no threshold. By default, we use \code{missing.rate = 0.1} to filter out variants with missing rate greater than 10\%.
#' @param ss.cutoff the minimum sample size (300 by default) to apply the MAF filter. This sample size is the sum of study samples and the benchmark samples of the same population as the study cohort.
#' @param maf to use the SNPs with ">= \code{maf}" if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise NaN is used by default for no MAF threshold.
#' @param hwe to use the SNPs with Hardy-Weinberg equilibrium p >= \code{hwe} if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise no hwe threshold. The default is 1e-6.
#' @param ... Arguments to be passed to other methods. 
#' @keywords PCA
#' @return a data frame with sample name, reported population, data resource (benchmark vs study cohort), the first four eigenvectors and the predicted population. 
#' @details
#' Using LD-pruned autosomal variants (by default), we calculate the eigenvectors and eigenvalues for principle component analysis (PCA). We use the benchmark samples as training dataset, and predict the population group for each sample in the study cohort based on the top four eigenvectors. Samples with discordant predicted and self-reported population groups are considered problematic. The function \code{PCACheck} performs the PCA analysis and identifies population outliers in study cohort.  
#' @importFrom stats predict
#' @export
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- PCACheck(seqfile, remove.samples=NULL, LDprune=TRUE, missing.rate=0.1)
#' res.pca <- QCresult(seqfile)$PCA
#' tail(res.pca)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

PCACheck <- function(seqfile, remove.samples = NULL, LDprune = TRUE, missing.rate = 0.1, ss.cutoff = 300, maf = 0.01, hwe = 1e-6, ...){
    
    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }

    message("calculating sample principle components ...")
    
    gfile <- SeqOpen(seqfile, readonly=TRUE)

    nds <- c("sample.id", "sample.annot", "snp.id") 
    allnds <- lapply(nds, function(x) read.gdsn(index.gdsn(gfile, x)))
    names(allnds) <- c("samples", "sampleanno", "snp.id")
    
    ## samples <- read.gdsn(index.gdsn(gfile, "sample.id"))
    ## snp.id <- read.gdsn(index.gdsn(gfile, "snp.id"))
    ## sampleanno <- read.gdsn(index.gdsn(gfile, "sample.annot"))

    ## remove samples for family related samples from 1kg, and other problematic samples from previous QC steps.
    study.pop <- unique(allnds$sampleanno[allnds$sampleanno$group == "study", "population"])
    if(length(study.pop) > 1) stop("Study samples should be single population, please prepare input file accordingly.")
    if(study.pop == "ASN") study.pop <- c("EAS", "SAS", "ASN")
    
    sample.relate <- allnds$sampleanno[allnds$sampleanno[,5]=="fam", 1]
    
    if(!is.null(remove.samples)){
        flag <- (allnds$sampleanno$group != "study" | allnds$sampleanno$population %in% study.pop) & ! allnds$samples %in% c(sample.relate, remove.samples)
    }else{
        flag <- (allnds$sampleanno$group != "study" | allnds$sampleanno$population %in% study.pop) & ! allnds$samples %in% sample.relate
    }
    sample.pca <- allnds$samples[flag]
    
    ## add hwe filter for variants when sample size >= 300.
    if (length(sample.pca) >= ss.cutoff){
        snp.hwe <- snpgdsHWE(gfile, sample.id=sample.pca)
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
    if (length(sample.pca) >= ss.cutoff){
        pca <- snpgdsPCA(gfile, sample.id=sample.pca, snp.id=allnds$snp.id[snp.idx], missing.rate=missing.rate, maf=maf, ...)
    }else{
        pca <- snpgdsPCA(gfile, sample.id=sample.pca, snp.id=allnds$snp.id[snp.idx], missing.rate=missing.rate, maf=NaN, ...)
    }   
    
    ## Make a data.frame and save result.
    res.pca <- data.frame(
        sample = pca$sample.id,
        pop = factor(allnds$sampleanno[,2])[match(pca$sample.id, allnds$sampleanno[,1])],
        type = factor(allnds$sampleanno[,5])[match(pca$sample.id, allnds$sampleanno[,1])],
        eval = pca$eigenval,
        ## pca$eigenvect,
        EV1 = pca$eigenvect[,1],    # the first eigenvector
        EV2 = pca$eigenvect[,2],    # the second eigenvector
        EV3 = pca$eigenvect[,3],
        EV4 = pca$eigenvect[,4],
        EV5 = pca$eigenvect[,5],
        EV6 = pca$eigenvect[,6],
        EV7 = pca$eigenvect[,7],
        EV8 = pca$eigenvect[,8],
        EV9 = pca$eigenvect[,9],
        EV10 = pca$eigenvect[,10], 
        stringsAsFactors = FALSE
    )
    ## Fst estimation
    ## pop.code <- factor(allnds$sampleanno[,2][match(pca$sample.id, allnds$sampleanno[,1])])
    ## a <- snpgdsFst(gfile, sample.id=sample.pca, population=pop.code, method="W&C84")
    
    ## ++++++++++++
    ## prediction 
    ## ++++++++++++
    
    ## use "svm" for classification, use the first 4 EigenVectors for prediction.
    ind.pop <- res.pca$type=="pop"
    model <- svm(pop ~ EV1 + EV2 + EV3 + EV4, data=res.pca[ind.pop, ], probability=FALSE, kernel="linear")
    pred.pop <- predict(model, res.pca[, 5:8], probability=FALSE)
    ## }
    res.pca <- cbind(res.pca[, c(1:3, 5:8)], pred.pop)

    ## for "ASN" samples, prediction of either "EAS" or "SAS" would be correct and summarized as "ASN".
    if("ASN" %in% levels(res.pca$pop)){
        res.pca[res.pca$pop == "ASN" & res.pca$pred.pop %in% c("EAS", "SAS"), "pred.pop"] <- "ASN"
    }
    ##res.pca <- cbind(res.pca, pred.pop, attributes(pred.pop)$probabilities)
    closefn.gds(gfile)

    ## return the SeqSQC file with updated gds file and QC result.
    a <- QCresult(seqfile)
    a$PCA <- res.pca
    outfile <- SeqSQC(gdsfile = gdsfile(seqfile), QCresult = a)
    return(outfile)
}
