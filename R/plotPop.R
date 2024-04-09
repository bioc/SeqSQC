plotPop <- function(res.qc, pc1="EV1", pc2="EV2", interactive=FALSE, pairedScatter=FALSE){
    if(pairedScatter){
        p <- ggpairs(
            res.qc,
            columns=c("EV1", "EV2", "EV3", "EV4"),
            aes(colour=pop, shape=type),
            upper="blank"
        )
    }else{
        cols <- match(c(pc1, pc2), colnames(res.qc))
        res.qc2 <- res.qc[, c(1:3, cols, 8)]
        colnames(res.qc2)[4:5] <- c("PC1", "PC2")
        p <- ggplot(data=res.qc2, aes(x=PC1, y=PC2))
        p <- p + geom_point(aes(colour=pop, shape=type))
        prob.pca <- res.qc2[res.qc2$pop != res.qc2$pred.pop, ]
        if(nrow(prob.pca) > 0){
            p <- p +
                geom_point(
                    data=prob.pca,
                    aes(x=PC1, y=PC2), shape=1, colour="red", size=4
                )
        }
        p <- p + theme_classic()
        ## p <- p + theme(legend.position = c(0.7, 0.8))
        p <- p + labs(x="Principle component 1",
                      y="Principle component 2",
                      col="Population",
                      shape="Type")
    }

    if(interactive){
        return(ggplotly(p))
    }else{
        return(p)
    }
}
