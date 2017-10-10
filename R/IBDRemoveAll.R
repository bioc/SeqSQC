IBDRemoveAll <- function(seqfile){

    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }
    if(!"IBD" %in% names(QCresult(seqfile))) stop("no IBD result.")

    sampleanno <- QCresult(seqfile)$sample.annot
    samples <- sampleanno$sample
    studyid <- sampleanno[sampleanno[,5] == "study", 1]
    res.ibd <- QCresult(seqfile)$IBD
    
    
    ## IBD results for all samples (including benchmark samples)
    res.study <- res.ibd
    ## res.study <- res.ibd[res.ibd$id1 %in% studyid & res.ibd$id2 %in% studyid,]

    ## check if any study sample is related with multiple samples
    ibd.related <- res.study[as.character(res.study$label) != as.character(res.study$pred.label) & res.study$kin > 0.08, ]
    ibd.pairs <- ibd.related[, 1:5]
    
    if(nrow(ibd.related) > 0){
        remove.multi <- c()
        repeat{
            ibd.multi <- table(c(as.character(ibd.related$id1), as.character(ibd.related$id2)))
            ibd.multi <- sort(ibd.multi, decreasing=TRUE)[1]
            if(ibd.multi==1){break}
            remove.multi <- rbind(remove.multi, c(names(ibd.multi), unname(ibd.multi)))
            ibd.related <- ibd.related[!ibd.related$id1 %in% remove.multi & !ibd.related$id2 %in% remove.multi, ]
        }
        ## remove the sample in the related pairs with higher missing rate. 
        if("MissingRate" %in% names(QCresult(seqfile))){
            mr <- QCresult(seqfile)$MissingRate
            mr <- mr[match(c(ibd.related$id1, ibd.related$id2), mr$sample), 1:2]
        }else{
            gfile <- SeqOpen(seqfile, readonly=TRUE)
            gt <- readex.gdsn(index.gdsn(gfile, "genotype"), list(NULL, samples %in% c(ibd.related$id1, ibd.related$id2)))
            mr <- apply(gt, 2, function(x) sum(! x %in% c(0, 1, 2))/length(x))
            mr <- data.frame(sample = c(ibd.related$id1, ibd.related$id2), missingRate = mr, stringsAsFactors=FALSE)
            closefn.gds(gfile)
        }    

        mr.rm <- function(x){
            mr.pair <- mr[match(x, mr$sample), 2]
            remove <- x[order(mr.pair, decreasing=TRUE)][1]
            return(remove)
        }
        remove.single <- unname(apply(ibd.related[, 1:2], 1, mr.rm))
        remove.ibd <- c(remove.multi[,1], remove.single)
    }else{
        remove.ibd=NULL
    }
    return(list(ibd.pairs=ibd.pairs, ibd.remove=remove.ibd))
}

        
