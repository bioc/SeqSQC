#' Plot the sample missing rate.
#' 
#' Plot sample missing rates calculated for each population.
#' @param seqfile SeqSQCclass input file with calculated missing rates.
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user perform sample QC using an rmarkdown script and output plot to html format.  
#' @import GGally
#' @import ggplot2
#' @import rbokeh
#' @export
#' @return the ggplot or interactive plot (if output is in html format) for sample missing rates.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gdsfile(seqfile) <- system.file("extdata", "example.gds", package="SeqSQC")
#' p <- plotMissingRate(seqfile, interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@@buffalo.edu}


plotMissingRate <- function(seqfile, interactive=FALSE){

    ## check
    if (!inherits(seqfile, "SeqSQCclass")){
        return("object should inherit from 'SeqSQCclass'.")
    }
    if(!"MissingRate" %in% names(QCresult(seqfile))) stop("no missing rate result.")

    sampleanno <- QCresult(seqfile)$sample.annot
    res.mr <- QCresult(seqfile)$MissingRate

    type <- sampleanno[match(res.mr$sample, sampleanno[,1]), 5]
    type[type=="fam"] <- "pop"
    res.mr$type <- factor(type, levels=c("pop", "study"), labels=c("benchmark", "study"))
    pop <- sampleanno[match(res.mr$sample, sampleanno[,1]), 2]
    res.mr$population <- factor(pop)
    rm(type, pop)
    
    ## by default, "study.pop" is a single population. If not single population, give error message.
    study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
    if(length(study.pop)>1) stop("Study samples should be single population.")
    
    ## only calculate "study.pop-specific" mean and sd, order results with pops. add all benchmark data in plotting.
    if(study.pop == "ASN"){
        res.pop <- res.mr[res.mr$population %in% c("EAS", "SAS", "ASN"), ]
    }else{
        res.pop <- res.mr[res.mr$population == study.pop, ]
    }
    
    b <- res.pop[res.pop$type == "benchmark", ]
    s <- res.pop[res.pop$type == "study", ]
    res.pop.ord <- rbind(b[order(b$population), ], s)
    res.pop.ord$sample <- factor(res.pop.ord$sample, levels=res.pop.ord$sample)
    res.pop.ord$outlier <- ifelse(res.pop.ord$missingRate > 0.1, "Yes", "No")
    res.pop.ord$problematic <- factor(res.pop.ord$outlier, levels=c("Yes", "No"))
    
    if(interactive){
        p1 <- suppressWarnings(figure(ylim=c(-0.1, 1))) %>% ly_points(missingRate, data=res.pop.ord, color = problematic, glyph = type, hover=list(sample, missingRate)) %>% ly_abline(h=0.1, type=2)
        return(p1)
    }else{
        p <- ggplot(data=res.pop.ord, aes(x=sample, y=missingRate))
        ##p <- p + scale_colour_manual(name="problematic",  values = c("Yes"="red", "No"="green"), drop=FALSE)
        p <- p + geom_point(aes(colour=problematic, shape=type))
        p <- p + scale_color_discrete(drop = FALSE)
        p <- p + geom_hline(yintercept=0.1, colour="grey")
        p <- p + theme(axis.text.x=element_blank())
        return(p)
    }
}
