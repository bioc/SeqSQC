#' Sample relationship check with SeqSQC object input file. 
#'
#' Function to calculate the IBD coefficients for all sample pairs and to predict related sample pairs in study cohort.
#'
#' @param seqfile SeqSQC object, which includes the merged gds file for study cohort and benchmark.
#' @param remove.samples a vector of sample names for removal from IBD calculation. Could be problematic samples identified from previous QC steps, or user-defined samples.
#' @param LDprune whether to use LD-pruned snp set. The default is TRUE.
#' @param kin.filter whether to use "kinship coefficient >= 0.08" as the additional criteria for related samples. The default is TRUE.
#' @param missing.rate to use the SNPs with "<= \code{missing.rate}" only; if NaN, no threshold. By default, we use \code{missing.rate = 0.1} to filter out variants with missing rate greater than 10\%.
#' @param ss.cutoff the minimum sample size (300 by default) to apply the MAF filter. This sample size is the sum of study samples and the benchmark samples of the same population as the study cohort.
#' @param maf to use the SNPs with ">= \code{maf}" if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise NaN is used by default for no MAF threshold.
#' @param hwe to use the SNPs with Hardy-Weinberg equilibrium p >= \code{hwe} if sample size defined in above argument is greater than \code{ss.cutoff}; otherwise no hwe threshold. The default is 1e-6.
#' @param ... Arguments to be passed to other methods. 
#' @keywords IBD
#' @return a data frame with sample names, the descent coefficients of k0, k1 and kinship, self-reported relationship and predicted relationship for each pair of samples.
#' @details Using LD-pruned variants (by default), we calculate the IBD coefficients for all sample pairs, and then predict related sample pairs in study cohort using the support vector machine (SVM) method with linear kernel and the known relatedness embedded in benchmark data as training set. \cr
#' Sample pairs with discordant self-reported and predicted relationship are considered as problematic. All predicted related pairs are also required to have coefficient of kinship >= 0.08 by default. The sample with higher missing rate in each related pair is selected for removal from further analysis by function of \code{IBDRemove}.
#' @import RColorBrewer
#' @import e1071
#' @import reshape2
#' @export
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- IBDCheck(seqfile, remove.samples=NULL, LDprune=TRUE, missing.rate=0.1)
#' res.ibd <- QCresult(seqfile)$IBD
#' tail(res.ibd)
#' @author Qian Liu \email{qliu7@@buffalo.edu}


IBDCheck <- function(seqfile, remove.samples = NULL, LDprune = TRUE, kin.filter = TRUE, missing.rate = 0.1, ss.cutoff = 300, maf = 0.01, hwe = 1e-6, ...){

    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }

    message("calculating pairwise IBD ...")
    
    gfile <- SeqOpen(seqfile, readonly=TRUE)
    samples <- read.gdsn(index.gdsn(gfile, "sample.id"))       
    sampleanno <- read.gdsn(index.gdsn(gfile, "sample.annot"))
    studyid <- sampleanno[sampleanno[,5] == "study", 1]
    snp.id <- read.gdsn(index.gdsn(gfile, "snp.id"))    

    ## sample filters. (remove prespecified "remove.samples", and only keep samples within the study population and benchmark data)
    study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
    if(length(study.pop) > 1) stop("Study samples should be single population, please prepare input file accordingly.")
    
    if(!is.null(remove.samples)){
        flag <- (sampleanno$group != "study" | sampleanno$population == study.pop) & !samples %in% remove.samples
    }else{
        flag <- sampleanno$group != "study" | sampleanno$population == study.pop
    }
    sample.ibd <- samples[flag]
    
    ## add hwe filter for variants when sample size >= 300.
    if (length(sample.ibd) >= ss.cutoff){
        snp.hwe <- snpgdsHWE(gfile, sample.id=sample.ibd)
        hwe.idx <- snp.hwe > hwe & !is.na(snp.hwe)
    }else{
        hwe.idx <- rep(TRUE, length(snp.id))
    }
    
    ## use LDpruned SNPs if LDprune == TRUE.
    if(LDprune){
        ld <- read.gdsn(index.gdsn(gfile, "snp.annot/LDprune"))
    }else{
        ld <- rep(TRUE, length(snp.id))
    }
    
    ## SNP filters in together. (HWE+LDprune)
    snp.idx <- ld & hwe.idx
    
    ## use maf filter if sample size >= 300.
    if (length(sample.ibd) >= ss.cutoff){ 
        IBD.res <- snpgdsIBDMoM(gfile, sample.id=sample.ibd, snp.id=snp.id[snp.idx], kinship=TRUE, maf=maf, missing.rate=missing.rate, ...)
    }else{
        IBD.res <- snpgdsIBDMoM(gfile, sample.id=sample.ibd, snp.id=snp.id[snp.idx], kinship=TRUE, maf=NaN, missing.rate=missing.rate, ...)
    }
    
    k0 <- IBD.res$k0
    k1 <- IBD.res$k1
    kin <- IBD.res$kinship
    
    ## We generate a big dataset for pairwise IBD.res result.
    sample.id.ibd <- IBD.res$sample.id
    rownames(kin) <- rownames(k1) <- rownames(k0) <- sample.id.ibd
    colnames(kin) <- colnames(k1) <- colnames(k0) <- sample.id.ibd
    idx <- melt(upper.tri(k0))[,3]
    k0m <- melt(k0)[idx,]
    k1m <- melt(k1)[idx,]
    kinm <- melt(kin)[idx,]
    res.ibd <- cbind(k0m, k1m[,3], kinm[,3])
    colnames(res.ibd) <- c("id1", "id2", "k0", "k1", "kin")        
    res.ibd$id1 <- as.character(res.ibd$id1)
    res.ibd$id2 <- as.character(res.ibd$id2)
    
    ## 1000g benchmark sample relation
    PO <- rbind(c("NA19238", "NA19240"),
                c("NA19239", "NA19240"),
                c("NA20868", "NA20871"),
                c("NA20886", "NA20898"))
    FS <- rbind(c("HG00581", "HG00635"),
                c("NA19713", "NA19985"))
    HF <- rbind(c("HG00119", "HG00124"),
                c("HG02353", "HG02363")
                )
    
    res.ibd$label <- "UN"
    res.ibd$label[paste(res.ibd$id1, res.ibd$id2) %in% paste(PO[,1], PO[,2])] <- "PO"
    res.ibd$label[paste(res.ibd$id1, res.ibd$id2) %in% paste(FS[,1], FS[,2])] <- "FS"
    res.ibd$label[paste(res.ibd$id1, res.ibd$id2) %in% paste(HF[,1], HF[,2])] <- "HF"
    res.ibd$label <- factor(res.ibd$label)
    
    ## index for 1000g benchmark data
    idx <- !(res.ibd$id1 %in% studyid | res.ibd$id2 %in% studyid)
    
    ## define relation centers
    centers <- data.frame(k0 = c(0, 0.25, 0.5, 1),
                          k1 = c(1, 0.5, 0.5, 0),
                          label = c("PO", "FS", "HF", "UN"),
                          stringsAsFactors=FALSE)
    
    ## if(method.pred=="svm"){
    ## use 1000g benchmark related samples and centers as training data. 
    ## svm: descriminate PO, FS, HF, FC, UN.
    res.svm <- rbind(res.ibd[idx, c("k0", "k1", "label")], centers[, c("k0", "k1", "label")])
    res.svm$label <- factor(res.svm$label)
    model <- svm(label ~ k0 + k1, data=res.svm, probability = FALSE, kernel="linear")
    ## predict the study sample relations, with probabilities.
    pred.svm <- predict(model, res.ibd[, c("k0", "k1")], probability = FALSE)
    pred.svm <- factor(pred.svm, levels=c("DU", levels(pred.svm)))
    pred.svm[res.ibd$kin > 0.46] <- "DU"
    ## use kin > 0.08 as additional filter.
    if(kin.filter){
        pred.svm[pred.svm != "UN" & res.ibd$kin < 0.08] <- "UN"
    }
    res.ibd <- data.frame(res.ibd[, 1:6], pred.label = pred.svm, stringsAsFactors=FALSE)
    res.ibd$pred.label <- as.character(res.ibd$pred.label)
    res.ibd$pred.label[res.ibd$pred.label != "UN"] <- "Related"
    res.ibd$pred.label <- factor(res.ibd$pred.label)
    ## }

    ## fully executable, could resume when updating.
    ## if(method.pred == "vglm"){
    ##     ## library(VGAM) ## importFrom VGAM vglm
    ##     ## use 1000g benchmark related samples and centers as training data.
    ##     ## vglm: descriminate PO, FS, HF, FC, UN.
    ##     res.vglm <- rbind(res.ibd[idx, c("k0", "k1", "label")], centers[, c("k0", "k1", "label")])
    ##     res.vglm$label <- factor(res.vglm$label)
    ##     model <- vglm(label ~ k0 + k1, data=res.vglm, family = multinomial)
    ##     probabilities <- predict(model, res.ibd[, c("k0", "k1")], type="response")
    ##     pred.vglm <- apply(probabilities, 1, which.max)
    ##     pred.vglm[which(pred.vglm=="1")] <- levels(res.vglm$label)[1]
    ##     pred.vglm[which(pred.vglm=="2")] <- levels(res.vglm$label)[2]
    ##     pred.vglm[which(pred.vglm=="3")] <- levels(res.vglm$label)[3]
    ##     pred.vglm[which(pred.vglm=="4")] <- levels(res.vglm$label)[4]

    ##     ## identify duplicates
    ##     pred.vglm[res.ibd$kin > 0.46] <- "DU"
    ##     pred.vglm <- factor(pred.vglm)
    ##     ## use kin > 0.08 as additional filter.
    ##     if(kin.filter){
    ##         pred.vglm[pred.vglm != "UN" & res.ibd$kin < 0.08] <- "UN"
    ##     }
    ##     res.ibd <- data.frame(res.ibd[, 1:6], pred.label = pred.vglm)
    ##     res.ibd$pred.label <- as.character(res.ibd$pred.label)
    ##     res.ibd$pred.label[res.ibd$pred.label != "UN"] <- "Related"
    ##     res.ibd$pred.label <- factor(res.ibd$pred.label)
    ## }
    closefn.gds(gfile)
    
    ## return the SeqSQC object with updated QC results.
    a <- QCresult(seqfile)
    a$IBD <- res.ibd
    outfile <- SeqSQC(gdsfile = gdsfile(seqfile), QCresult = a)
    return(outfile)


}



