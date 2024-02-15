plotInbreeding <- function(res.qc, interactive = FALSE, sdcoef = 5){
    res.pop.ord <- res.qc
    inb <- res.pop.ord$inbreeding
    m <- mean(inb)
    sd <- sd(inb)
    sd1 <- m + sdcoef*sd
    sd2 <- m - sdcoef*sd
    res.pop.ord$outlier <- ifelse(res.pop.ord$inbreeding > sd1 | res.pop.ord$inbreeding < sd2,
                                  "Yes", "No")
    res.pop.ord$problematic <- factor(res.pop.ord$outlier, levels=c("Yes", "No"))
    
    p <- ggplot(data=res.pop.ord, aes(x=sample, y=inbreeding))
    p <- p + geom_point(aes(shape=type, colour=problematic))
    p <- p + scale_colour_discrete(drop = FALSE)
    p <- p + geom_hline(yintercept=c(sd1, sd2), colour="grey")
    p <- p + theme_classic()
    p <- p + labs(col="Problematic",
                  shape="Type",
                  x="Sample",
                  y="Autosomal inbreeding coefficient")
    ## p <- p + theme(legend.position = c(0.1, 0.2))
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())  

    if(interactive){
        return(ggplotly(p))
    }else{
        return(p)
    }
}


