#' Plot the QC results for specific QC steps. 
#' 
#' Plot QC results. 
#' @param seqfile SeqSQC object with QC results.
#' @param QCstep which QC step the user want to do plotting. 
#' @param interactive whether to generate interactive plot. Recommend to use \code{interactive = TRUE} if user perform sample QC using an rmarkdown script and output plot to html format.
#' @param sdcoef for inbreeding outlier check, how many standard deviation we need for identification of inbreeding outliers. The default is 5.
#' @param pc1 the eigenvector on x axis for PCA result. The default is "EV1" for eigenvector 1.
#' @param pc2 the eigenvector on y axis for PCA result. The default is "EV2" for eigenvector 2. 
#' @param pairedScatter for PCA result, whether to plot the paired scatterplot for the first 4 PC axes.
#' @param ... Arguments to be passed to other methods. 
#' @import ggplot2
#' @import RColorBrewer
#' @import GGally
#' @import rbokeh
#' @export
#' @return the ggplot or interactive plot (if output is in html format) for specific QC result.
#' @examples
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' p <- plotQC(seqfile, QCstep="PCA", interactive=FALSE)
#' p
#' @author Qian Liu \email{qliu7@buffalo.edu}

plotQC <- function(seqfile, QCstep=NULL, interactive=FALSE, sdcoef=5, pc1="EV1", pc2="EV2", pairedScatter=FALSE, ...){

    res.qc <- plotQCData(seqfile, QCstep = QCstep)
    
    ## plotting -- MissingRate
    if(QCstep == "MissingRate"){
        p <- plotMissingRate(res.qc, interactive = interactive)
    }else if(QCstep == "SexCheck"){
        p <- plotSexCheck(res.qc, interactive = interactive)
    }else if(QCstep == "Inbreeding"){
        p <- plotInbreeding(res.qc, interactive = interactive, sdcoef = sdcoef)
    }else if(QCstep == "IBD"){
        p <- plotIBD(res.qc, interactive = interactive)
    }else if(QCstep == "PCA"){
     p <- plotPop(res.qc, interactive = interactive, ...)
    }
    return(p)
}
