plotQCData <- function(seqfile, QCstep=NULL){

    ## input check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }
    if(!QCstep %in% names(QCresult(seqfile))) stop(paste("no QC result for", QCstep))

    if(QCstep == "MissingRate"){
        res.qc <- QCresult(seqfile)$MissingRate
    }else if(QCstep == "SexCheck"){
        res.qc <- QCresult(seqfile)$SexCheck
    }else if(QCstep == "Inbreeding"){
        res.qc <- QCresult(seqfile)$Inbreeding
    }else if(QCstep == "IBD"){
        res.qc <- QCresult(seqfile)$IBD
    }else if(QCstep == "PCA"){
        res.qc <- QCresult(seqfile)$PCA
    }
    sampleanno <- QCresult(seqfile)$sample.annot
    
    if(QCstep %in% c("MissingRate", "SexCheck", "Inbreeding")){
        type <- sampleanno[match(res.qc$sample, sampleanno[,1]), 5]
        type[type=="fam"] <- "pop"
        res.qc$type <- factor(type, levels=c("pop", "study"), labels=c("benchmark", "study"))
        pop <- sampleanno[match(res.qc$sample, sampleanno[,1]), 2]
        res.qc$population <- factor(pop)
        rm(type, pop)
        
        ## by default, "study.pop" is a single population for missing
        ## rate, sexcheck and inbreeding. If not single population,
        ## give error message.
        study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
        if(length(study.pop)>1) stop("Study samples should be single population.")
        
        ## population specific calculation for missingrate, sexcheck and inbreeding. 
        if(study.pop == "ASN"){
            res.pop <- res.qc[res.qc$population %in% c("EAS", "SAS", "ASN"), ]
        }else{
            res.pop <- res.qc[res.qc$population == study.pop, ]
        }

        b <- res.pop[res.pop$type == "benchmark", ]
        s <- res.pop[res.pop$type == "study" & res.pop$population == study.pop, ]
        res.pop.ord <- rbind(b[order(b$population), ], s)
        res.pop.ord$sample <- factor(res.pop.ord$sample, levels=res.pop.ord$sample)
        res.qc <- res.pop.ord
        
    }else if(QCstep == "IBD"){
        studyid <- sampleanno[sampleanno[,5] == "study", 1]
        
        res.qc$label <- as.character(res.qc$label)
        res.qc$relation <- factor(res.qc$label, levels=c("DU", "PO", "FS", "HF", "UN"))
        
        ind.study <- res.qc$id1 %in% studyid | res.qc$id2 %in% studyid
        names(res.qc)[6] <- "type"
        res.qc$type[ind.study] <- "study"
        res.qc$type[!ind.study] <- "benchmark"
        
    }else if(QCstep == "PCA"){
        res.qc$type <- as.character(res.qc$type)
        res.qc$type[res.qc$type == "pop"] <- "benchmark"
        res.qc$type <- factor(res.qc$type)
    }
    res.qc
}
