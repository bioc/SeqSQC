#' Plot the PC axes for samples.
#' 
#' Plot principle components for all sample.
#' @param seqfile SeqSQC object with PCA results.
#' @param pc1 the eigenvector on x axis. The default is "EV1" for eigenvector 1.
#' @param pc2 the eigenvector on y axis. The default is "EV2" for eigenvector 2. 
#' @param pairedScatter whether to plot the paired scatterplot for the first 4 PC axes.
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user perform sample QC using an rmarkdown script and output plot to html format.  
#' @export
#' @return the ggplot or interactive plot (if output is in html format) for principle components.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' p <- plotPop(seqfile, interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@@buffalo.edu}

plotPop <- function(seqfile, pc1="EV1", pc2="EV2", interactive=FALSE, pairedScatter=FALSE){

    ## check
    if (!inherits(seqfile, "SeqSQC")){
        return("object should inherit from 'SeqSQC'.")
    }
    if(!"PCA" %in% names(QCresult(seqfile))) stop("no PCA result.")

    res.pca <- QCresult(seqfile)$PCA
    res.pca$type <- as.character(res.pca$type)
    res.pca$type[res.pca$type == "pop"] <- "benchmark"
    res.pca$type <- factor(res.pca$type)

    if(interactive){
        if(pairedScatter){
            tools <- c("pan", "wheel_zoom", "box_zoom", "box_select", "reset")
            nms <- expand.grid(colnames(res.pca)[5:8], rev(colnames(res.pca)[4:7]), stringsAsFactors = FALSE)
            splom.list <- vector("list", 16)
            for(i in seq_len(nrow(nms))) {
                splom.list[[i]] <- figure(width = 200, height = 200, tools = tools,
                                          xlab = nms$Var1[i], ylab = nms$Var2[i]) %>%
                    ly_points(nms$Var1[i], nms$Var2[i], data = res.pca,
                              color = pop, size = 5, legend = FALSE, hover=list(sample))
            }
            splom.list[[4]] <- figure(width = 200, height = 200,
                                      xlab = nms$Var1[i], ylab = nms$Var2[i]) %>%
                ly_points(nms$Var1[i], nms$Var2[i], data = res.pca,
                          color = pop, visible=FALSE)
            p <- grid_plot(splom.list, ncol = 4, same_axes = TRUE, link_data = TRUE)
            return(p)
        }else{
            cols <- match(c(pc1, pc2), colnames(res.pca))
            res.pca2 <- res.pca[, c(1:3, cols, 8)]
            colnames(res.pca2)[4:5] <- c("PC1", "PC2")
            p <- figure() %>% ly_points(PC1, PC2, data=res.pca2, color=pop, glyph=type, hover=list(sample, pop, pred.pop))
        }
        return(p)
    }else{
        if(pairedScatter){
            p <- ggpairs(res.pca, columns=c("EV1", "EV2", "EV3", "EV4"), aes(colour=pop, shape=type), upper="blank")
        }else{
            cols <- match(c(pc1, pc2), colnames(res.pca))
            res.pca2 <- res.pca[, c(1:3, cols, 8)]
            colnames(res.pca2)[4:5] <- c("PC1", "PC2")
            p <- ggplot(data=res.pca2, aes(x=PC1, y=PC2))
            p <- p + geom_point(aes(colour=pop, shape=type))
            prob.pca <- res.pca2[res.pca2$pop != res.pca2$pred.pop, ]
            if(nrow(prob.pca) > 0){
                p <- p + geom_point(data=prob.pca, aes(x=PC1, y=PC2), shape=1, colour="red", size=4)
            }
        }
        return(p)
    }
}
