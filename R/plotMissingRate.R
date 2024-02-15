plotMissingRate <- function(res.qc, interactive = FALSE){

    res.pop.ord <- res.qc
    res.pop.ord$outlier <- ifelse(res.pop.ord$missingRate > 0.1, "Yes", "No")
    res.pop.ord$problematic <- factor(res.pop.ord$outlier, levels=c("Yes", "No"))
    
    p <- ggplot(data=res.pop.ord, aes(x=sample, y=missingRate))
    p <- p + geom_point(aes(colour=problematic, shape=type))
    p <- p + scale_color_discrete(drop = FALSE)
    p <- p + geom_hline(yintercept=0.1, colour="grey")
    p <- p + theme(axis.text.x=element_blank())

    if(interactive){
        return(ggplotly(p))
    }else{
        return(p)
    }

}

