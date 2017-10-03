#' Plot the inbreeding coefficients.
#'
#' Plot inbreeding coefficients calculated from autosomal variants for each population.
#' @param seqfile SeqSQCclass input file with inbreeding result.
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user perform sample QC using an rmarkdown script and output plot to html format.
#' @param sdcoef how many standard deviation we need for identification of inbreeding outliers. The default is 5.
#' @export
#' @return @return the ggplot or interactive plot (if output is in html format) for inbreeding coefficients.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQCclass(gdsfile = gfile, QCresult = QCresult(seqfile))
#' p <- plotInbreeding(seqfile, interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@@buffalo.edu}

plotInbreeding <- function(seqfile, interactive=FALSE, sdcoef=5){

    ## check
    if (!inherits(seqfile, "SeqSQCclass")){
        return("object should inherit from 'SeqSQCclass'.")
    }
    if(!"Inbreeding" %in% names(QCresult(seqfile))) stop("No inbreeding result.")

    sampleanno <- QCresult(seqfile)$sample.annot
    res.inb <- QCresult(seqfile)$Inbreeding

    ## define "res.inb$type" as "study/benchmark".
    type <- sampleanno[match(res.inb$sample, sampleanno[,1]), 5]
    type[type=="fam"] <- "pop"
    res.inb$type <- factor(type, levels=c("pop", "study"), labels=c("benchmark", "study"))
    pop <- sampleanno[match(res.inb$sample, sampleanno[,1]), 2]
    res.inb$population <- factor(pop)
    rm(type, pop)
    
    ## by default, "study.pop" is a single population. If not single population, give error message.
    study.pop <- unique(sampleanno[sampleanno$group == "study", "population"])
    if(length(study.pop)>1) stop("Study samples should be single population.")
    
    ## only calculate "study.pop-specific" mean and sd, order results with pops. add all benchmark data in plotting.
    if(study.pop == "ASN"){
        res.pop <- res.inb[res.inb$population %in% c("EAS", "SAS", "ASN"), ]
    }else{
        res.pop <- res.inb[res.inb$population == study.pop, ]
    }
    
    inb <- res.pop$inbreeding
    m <- mean(inb)
    sd <- sd(inb)

    b <- res.pop[res.pop$type == "benchmark", ]
    s <- res.pop[res.pop$type == "study" & res.pop$population == study.pop, ]
    res.pop.ord <- rbind(b[order(b$population), ], s)
    res.pop.ord$sample <- factor(res.pop.ord$sample, levels=res.pop.ord$sample)

    sd1 <- m + sdcoef*sd
    sd2 <- m - sdcoef*sd
    res.pop.ord$outlier <- ifelse(res.pop.ord$inbreeding > sd1 | res.pop.ord$inbreeding < sd2, "Yes", "No")
    res.pop.ord$problematic <- factor(res.pop.ord$outlier, levels=c("Yes", "No"))
    ## iterative calculation of mean and sd. 
    ## repeat{
    ##     t <- which(inb>m+sdcoef*s | inb<m-sdcoef*s)
    ##     if (length(t)==0) break
    ##     if (mean(inb[-t])==m & sd(inb[-t])==s) break
    ##     m <- mean(inb[-t])
    ##     s <- sd(inb[-t])
    ## }

    if(interactive){
        s1 <- max(max(res.pop.ord$inbreeding), sd1) + 0.1
        s2 <- min(min(res.pop.ord$inbreeding), sd2) - 0.1
        p1 <- figure(ylim=c(s2, s1)) %>% ly_points(inbreeding, data=res.pop.ord, color = problematic, glyph = type, hover=list(sample, inbreeding)) %>% ly_abline(h=sd1, type=2) %>% ly_abline(h=sd2, type=2)
        return(p1)
    }else{
        p <- ggplot(data=res.pop.ord, aes(x=sample, y=inbreeding))
        p <- p + geom_point(aes(shape=type, colour=problematic))
        p <- p + scale_colour_discrete(drop = FALSE)
        p <- p + geom_hline(yintercept=c(sd1, sd2), colour="grey")
        p <- p + theme(axis.text.x=element_blank())
        return(p)
    }
}

