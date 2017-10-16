plotPop <- function(res.qc, pc1="EV1", pc2="EV2", interactive=FALSE, pairedScatter=FALSE){
    if(interactive){
        if(pairedScatter){
            tools <- c("pan", "wheel_zoom", "box_zoom", "box_select", "reset")
            nms <- expand.grid(
                colnames(res.qc)[5:8],
                rev(colnames(res.qc)[4:7]),
                stringsAsFactors = FALSE
            )
            splom.list <- vector("list", 16)
            for(i in seq_len(nrow(nms))) {
                splom.list[[i]] <- figure(
                    width = 200, height = 200, tools = tools,
                    xlab = nms$Var1[i], ylab = nms$Var2[i]) %>%
                    ly_points(
                        nms$Var1[i], nms$Var2[i], data = res.qc,
                        color = pop, size = 5, legend = FALSE, hover=list(sample))
            }
            splom.list[[4]] <- figure(
                width = 200, height = 200,
                xlab = nms$Var1[i], ylab = nms$Var2[i]) %>%
                ly_points(nms$Var1[i], nms$Var2[i], data = res.qc,
                          color = pop, visible=FALSE)
            p <- grid_plot(
                splom.list, ncol = 4, same_axes = TRUE, link_data = TRUE
            )
            return(p)
        }else{
            cols <- match(c(pc1, pc2), colnames(res.qc))
            res.qc2 <- res.qc[, c(1:3, cols, 8)]
            colnames(res.qc2)[4:5] <- c("PC1", "PC2")
            p <- figure() %>%
                    ly_points(
                        PC1, PC2, data=res.qc2, color=pop,
                        glyph=type, hover=list(sample, pop, pred.pop)
                    )
        }
        return(p)
    }else{
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
        }
        return(p)
    }
}
