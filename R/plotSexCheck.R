plotSexCheck <- function(res.qc, interactive=FALSE){
    if (any(res.qc$sex == "")) {
        warning("Please provide sex info for all study samples!")
        res.qc <- res.qc[res.qc$sex != "", ]
    }
    res.pop.ord <- res.qc
    prob.ind <- res.pop.ord$sex == "female" & res.pop.ord$pred.sex == "male" | res.pop.ord$sex == "male" & res.pop.ord$pred.sex == "female"
    prob.sex <- res.pop.ord[prob.ind, ]
    
    if(interactive){
        p1 <- figure(ylim=c(-2,2)) %>%
            ly_points(
                sexinb, data=res.pop.ord, color = sex, glyph = type,
                hover=list(sample, sex, pred.sex)
            ) %>%
            x_axis(label="Sample") %>%
            y_axis(label="Sex inbreeding coefficient") %>%
            ly_abline(h=c(0.2, 0.8), type=2)
        return(p1)
    }else{
        p <- ggplot(data=res.pop.ord, aes(x=sample, y=sexinb))
        p <- p + geom_point(aes(colour=sex, shape=type))
        p <- p + scale_colour_discrete(drop = FALSE)
        p <- p + labs(col="Self-reported gender",
                      shape="Type",
                      x="Sample",
                      y="Sex inbreeding coefficient")
        p <- p + geom_hline(yintercept=c(0.2, 0.8), colour="grey")
        if(nrow(prob.sex) > 0){
            p <- p + geom_point(data=prob.sex, aes(x=sample, y=sexinb), shape=1, colour="red", size=4)
        }
        p <- p + theme_classic()
        p <- p + theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())  
        ## p <- p + theme(legend.position = c(0.2, 0.5))
        return(p)
    }
}
