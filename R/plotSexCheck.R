#' Plot the X chromosome inbreeding coefficients.
#' 
#' Plot the inbreeding coefficients calculated from variants on chromosome X for each populaton.
#' @param seqfile SeqSQCclass input file with sexcheck results.
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user perform sample QC using an rmarkdown script and output plot to html format. 
#' @export
#' @return the ggplot or interactive plot (if output is in html format) for sex inbreeding coefficients.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gdsfile(seqfile) <- system.file("extdata", "example.gds", package="SeqSQC")
#' p <- plotSexCheck(seqfile, interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@@buffalo.edu}

plotSexCheck <- function(seqfile, interactive=FALSE){

    ## check
    if (!inherits(seqfile, "SeqSQCclass")){
        return("object should inherit from 'SeqSQCclass'.")
    }
    if(!"SexCheck" %in% names(QCresult(seqfile))) stop("no sexcheck result.")

    sampleanno <- QCresult(seqfile)$sample.annot
    res.sex <- QCresult(seqfile)$SexCheck

    type <- sampleanno[match(res.sex$sample, sampleanno[,1]), 5]
    type[type=="fam"] <- "pop"
    res.sex$type <- factor(type, levels=c("pop", "study"), labels=c("benchmark", "study"))
    pop <- sampleanno[match(res.sex$sample, sampleanno[,1]), 2]
    res.sex$population <- factor(pop)
    rm(type, pop)
    
    ## by default, "study.pop" is a single population. If not single population, give error message.
    study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
    if(length(study.pop)>1) stop("Study samples should be single population.")

    ## only calculate "study.pop-specific" mean and sd, order results with pops. add all benchmark data in plotting.
    if(study.pop == "ASN"){
        res.pop <- res.sex[res.sex$population %in% c("EAS", "SAS", "ASN"), ]
    }else{
        res.pop <- res.sex[res.sex$population == study.pop, ]
    }
    
    b <- res.pop[res.pop$type == "benchmark", ]
    s <- res.pop[res.pop$type == "study", ]
    res.pop.ord <- rbind(b[order(b$population), ], s)
    res.pop.ord$sample <- factor(res.pop.ord$sample, levels=res.pop.ord$sample)

    ## add dashed circle for samples with gender mismatch
    prob.ind <- res.pop.ord$sex == "female" & res.pop.ord$pred.sex == "male" | res.pop.ord$sex == "male" & res.pop.ord$pred.sex == "female"
    prob.sex <- res.pop.ord[prob.ind, ]
    
    if(interactive){
        p1 <- figure(ylim=c(-2,2)) %>% ly_points(sexinb, data=res.pop.ord, color = sex, glyph = type, hover=list(sample, sex, pred.sex)) %>% ly_abline(h=c(0.2, 0.8), type=2)
        return(p1)
    }else{
        p <- ggplot(data=res.pop.ord, aes(x=sample, y=sexinb))
        p <- p + geom_point(aes(colour=sex, shape=type))
        p <- p + scale_colour_discrete(drop = FALSE)
        p <- p + labs(col="self-reported gender")
        p <- p + geom_hline(yintercept=c(0.2, 0.8), colour="grey")
        if(nrow(prob.sex) > 0){
            p <- p + geom_point(data=prob.sex, aes(x=sample, y=sexinb), shape=1, colour="red", size=4)
        }
        p <- p + theme(axis.text.x=element_blank())  
        return(p)
    }
}
