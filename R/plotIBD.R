#' Plot the IBD coefficients for all sample pairs.
#' 
#' Plot IBD coefficients (k0 and k1) calculated from autosomal variants for all samples.
#' @param seqfile SeqSQCclass input file with IBD results.
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user performs sample QC using an rmarkdown script and output plots to html format. 
#' @export
#' @return the ggplot or interactive plot (if output is in html format) for IBD coefficients.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gdsfile(seqfile) <- system.file("extdata", "example.gds", package="SeqSQC")
#' p <- plotIBD(seqfile, interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@@buffalo.edu}

plotIBD <- function(seqfile, interactive=FALSE){

    ## check
    if (!inherits(seqfile, "SeqSQCclass")){
        return("object should inherit from 'SeqSQCclass'.")
    }
    if(!"IBD" %in% names(QCresult(seqfile))) stop("no IBD result.")

    sampleanno <- QCresult(seqfile)$sample.annot
    res.ibd <- QCresult(seqfile)$IBD
    studyid <- sampleanno[sampleanno[,5] == "study", 1]
    
    centers <- data.frame(k0=c(0, 0, 0.25, 0.5
                             ## , 0.75
                             , 1),
                          k1=c(0, 1, 0.5, 0.5
                               ## , 0.25
                             , 0),
                          label=c("DU", "PO", "FS", "HF"
                                  ## , "FC"
                                , "UN"))
    centers$label <- factor(centers$label, levels=centers$label)
    res.ibd$label <- as.character(res.ibd$label)
    res.ibd$relation <- factor(res.ibd$label, levels=as.character(centers$label))

    ind.study <- res.ibd$id1 %in% studyid | res.ibd$id2 %in% studyid
    res.ibd$label[ind.study] <- "study"
    res.ibd$label[!ind.study] <- "benchmark"
    names(res.ibd)[6] <- "type"

    prob.ind <- as.character(res.ibd$relation) == "UN" & as.character(res.ibd$pred.label) == "Related"
    prob.ibd <- res.ibd[prob.ind, ]
    
    if(interactive){
        p <- figure() %>% ly_points(k0, k1, data=res.ibd, hover=list(relation, pred.label, id1, id2), color=relation) %>% ly_points(k0, k1, data=centers, hover=list(label), color="grey", glyph=label, legend=FALSE)
        return(p)
    }else{
        p <- ggplot(data=res.ibd, aes(x=k0, y=k1))
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

