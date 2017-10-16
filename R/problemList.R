#' Generate the problematic sample list.
#'
#' generate the problematic sample list from QC steps that have been done, and provide each problematic sample with a reason for removal (high missing rate, gender mismatch, inbreeding outlier, cryptic relationship or population outlier).
#' @param seqfile SeqSQC object with sample QC results. 
#' @export
#' @return a data frame with 2 columns: \code{sample} for problematic sample name, and \code{remove.reason} for the reason of removing the sample. 
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' problemList(seqfile)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

problemList <- function(seqfile){

    sampleanno <- QCresult(seqfile)$sample.annot
    samples <- sampleanno$sample
    studyid <- sampleanno[sampleanno[,5] == "study", 1]

    if("MissingRate" %in% names(QCresult(seqfile))){
        res.mr <- QCresult(seqfile)$MissingRate
        res.study <- res.mr[res.mr$sample %in% studyid, ]
        prob.mr <- res.study[res.study$outlier == "Yes", ]
        remove.mr <- prob.mr[,1]
    }else{
        remove.mr <- NULL
    }
        
    if("SexCheck" %in% names(QCresult(seqfile))){
        res.sexcheck <- QCresult(seqfile)$SexCheck
        prob.sex <- res.sexcheck[res.sexcheck$sex == "female" & res.sexcheck$pred.sex == "male" | res.sexcheck$sex == "male" & res.sexcheck$pred.sex == "female", ]
        remove.sex <- prob.sex[,1]
    }else{
        remove.sex <- NULL
    }

    if("Inbreeding" %in% names(QCresult(seqfile))){
    res.inb <- QCresult(seqfile)$Inbreeding
    prob.inb <- res.inb[res.inb$outlier.5sd == "Yes" & !is.na(res.inb$outlier.5sd), ]
    remove.inb <- prob.inb[,1]
    }else{
        remove.inb <- NULL
    }

    if("IBD" %in% names(QCresult(seqfile))){
        res.ibd <- QCresult(seqfile)$IBD
        prob.ibd <- IBDRemove(seqfile)
        ibd.pairs <- prob.ibd$ibd.pairs
        if(nrow(ibd.pairs) > 0){
            ibd.pairs <- paste(ibd.pairs$id1, ibd.pairs$id2, sep=":")
        }else{
            ibd.pairs = NULL
        }
        remove.ibd <- prob.ibd$ibd.remove
    }else{
        ibd.pairs <- NULL
        remove.ibd <- NULL
    }
    
    if("PCA" %in% names(QCresult(seqfile))){
        res.pca <- QCresult(seqfile)$PCA
        res.study <- res.pca[res.pca$sample %in% studyid, ]
        remove.pca <- res.study[res.study$pop != res.study$pred.pop, 1]
    }else{
        remove.pca <- NULL
    }
    
    prob.list <- data.frame(
        sample=c(remove.mr, remove.sex, remove.inb, ibd.pairs, remove.pca),
        ## remove=c(remove.mr, remove.sex, remove.inb, remove.ibd, remove.pca),
        remove.reason=c(
            rep("high missing rate", length(remove.mr)),
            rep("gender mismatch", length(remove.sex)),
            rep("inbreeding outlier", length(remove.inb)),
            rep("cryptic relationship", length(ibd.pairs)),
            rep("population outlier", length(remove.pca))),
        stringsAsFactors=FALSE
    )
    remove.list <- data.frame(sample = c(remove.mr, remove.sex, remove.inb, remove.ibd, remove.pca))
    if(nrow(prob.list) != 0){
        a <- QCresult(seqfile)
        a$problem.list <- list(prob.list = prob.list, remove.list = remove.list)
    }
    return(a$problem.list)
}

