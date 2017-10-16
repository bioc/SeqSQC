plotIBD <- function(res.qc, interactive = FALSE){

    centers <- data.frame(k0=c(0, 0, 0.25, 0.5, 1),
                          k1=c(0, 1, 0.5, 0.5, 0),
                          label=factor(c("DU", "PO", "FS", "HF", "UN"))
                          )
    prob.ind <- as.character(res.qc$relation) == "UN" & as.character(res.qc$pred.label) == "Related"
    prob.ibd <- res.qc[prob.ind, ]

    if(interactive){
        p <- figure() %>% ly_points(k0, k1, data=res.qc, hover=list(relation, pred.label, id1, id2), color=relation) %>% ly_points(k0, k1, data=centers, hover=list(label), color="grey", glyph=label, legend=FALSE)
        return(p)
    }else{
        p <- ggplot(data=res.qc, aes(x=k0, y=k1))
        p <- p + geom_point(aes(shape=type, colour=relation))
        p <- p + scale_colour_discrete(drop = FALSE)
        p <- p + geom_point(aes(colour=label), shape=3, data=centers)
        cols <- c(brewer.pal(n=9, name="Set1")[c(1,3,5,8)], "black")
        p <- p + scale_color_manual(values=cols)
        if(nrow(prob.ibd) > 0){
            p <- p + geom_point(data=prob.ibd, aes(x=k0, y=k1), shape=1, colour="red", size=4)
        }
        return(p)
    }
}



