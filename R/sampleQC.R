#' The wrap-up function for sample QC of sequencing/GWAS data.
#' 
#' A wrap-up function for sample QC. It reads in the variant genotypes
#' in vcf/PLINK format, merges study cohort with benchmark data, and
#' performs sample QC for the merged dataset.
#' 
#' @param vfile vcf or PLINK input file (ped/map/bed/bim/fam with same
#'     basename). The default is NULL. Vfile could be a vector of
#'     character strings, see details. Could also take file in
#'     \code{SeqSQC} object generated from \code{LoadVfile}.
#' @param output a character string for name of merged data of SeqSQC
#'     object. The \code{dirname(output)} would be used as the
#'     directory to save the QC result and plots. The default is
#'     \code{sampleqc} in the working directory.
#' @param capture.region the BED file of sequencing capture
#'     regions. The default is NULL. For exome-sequencing data, the
#'     capture region file must be provided.
#' @param sample.annot sample annotation file with 3 columns (with
#'     header) in the order of sample id, sample population and sex
#'     info. The default is NULL.
#' @param LDprune whether to use LD-pruned snp set. The default is
#'     TRUE.
#' @param vfile.restrict whether the input vcf or plink file has
#'     already been restricted by capture region. The default is
#'     FALSE.
#' @param slide.max.bp the window size of SNPs when calculating
#'     linkage disequilibrium. The default is 5e+05.
#' @param ld.threshold the r^2 threshold for LD-based SNP pruning if
#'     \code{LDprune = TRUE}. The default is 0.3.
#' @param format.data the data source. The default is \code{NGS} for
#'     sequencing data.
#' @param format.file the data format. The default is \code{vcf}.
#' @param QCreport Whether to generate the sample QC report in html
#'     format.
#' @param out.report the file name for the sample QC report. The
#'     default is \code{report.html}.
#' @param interactive whether to generate interactive plots in the
#'     sample QC report if \code{QCreport = TRUE}.
#' @param results whether to write out the results for each QC steps
#'     in .txt files. The default is TRUE.
#' @param plotting whether to output the plots for each QC steps in
#'     .pdf files. the default is TRUE.
#' @param ... Arguments to be passed to other methods.
#' @export
#' @return a SeqSQC object with the filepath to the gds file which
#'     stores the genotype, the summary of samples and variants, and
#'     the QCresults including the sample annotation information and
#'     all QC results.

#' @details For \code{vfile} with more than one file names,
#'     \code{sampleQC} will merge all dataset together if they all
#'     contain the same samples. It is useful to combine
#'     genetic/genomic data together if VCF data is divided by
#'     chromosomes. \cr There are 3 columns in \code{sample.annot}
#'     file. col 1 is \code{sample} with sample ids, col 2 is
#'     \code{population} with values of "AFR/EUR/ASN/EAS/SAS", col 3
#'     is \code{gender} with values of "male/female".
#' @importFrom utils write.table
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "example_sub.vcf", package="SeqSQC")
#' sample.annot <- system.file("extdata", "sampleAnnotation.txt", package="SeqSQC")
#' cr <- system.file("extdata", "CCDS.Hs37.3.reduced_chr1.bed", package="SeqSQC")
#' outfile <- file.path(tempdir(), "testWrapUp")
#' seqfile <- sampleQC(vfile = infile, output = outfile, capture.region = cr,
#' sample.annot = sample.annot, format.data = "NGS", format.file = "vcf",
#' QCreport = TRUE, out.report="report.html", interactive = TRUE)
#' ## save(seqfile, file="seqfile.RData")
#'
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' seqfile <- sampleQC(sfile = seqfile, output = outfile, QCreport = FALSE,
#' out.report="report.html", interactive = TRUE)
#' }
#' @author Qian Liu \email{qliu7@buffalo.edu}

sampleQC <- function(vfile = NULL, output="sampleqc",
                     capture.region = NULL, sample.annot = NULL,
                     LDprune = TRUE, vfile.restrict = FALSE,
                     slide.max.bp = 5e+05, ld.threshold = 0.3,
                     format.data = "NGS", format.file = "vcf",
                     QCreport = TRUE, out.report="report.html",
                     interactive = TRUE, results=TRUE, plotting=TRUE,
                     ...){

    ## check input
    if(inherits(vfile, "SeqSQC")){
        seqfile <- vfile
        format.data <- NULL
        format.file <- NULL
        vfile.restrict <- NULL
    }else{
        seqfile <- LoadVfile(vfile = vfile, output = output, capture.region = capture.region, sample.annot = sample.annot, ...)
    }
    fn <- gdsfile(seqfile)
    print(paste("gds file generated:", fn))
  
    sampleanno <- QCresult(seqfile)$sample.annot
    samples <- sampleanno$sample
    studyid <- sampleanno[sampleanno[,5] == "study", 1]
    
    ##############
    ## sample missing rate
    
    seqfile <- MissingRate(seqfile)
    res.mr <- QCresult(seqfile)$MissingRate
    res.study <- res.mr[res.mr$sample %in% studyid, ]
    prob.mr <- res.study[res.study$outlier == "Yes", ]
    remove.mr <- prob.mr[,1]
    remove.samples <- unique(remove.mr)
    if(results){
        write.table(res.mr, file=paste0(dirname(output), "/result.missingrate.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    }
    if(plotting){
        p <- plotQC(seqfile, "MissingRate")
        ggsave(filename = paste0(dirname(output), "/plot.missingrate.pdf"), p)
    }
    
    ##############
    ## Sex Check

    ## check chrX variants
    gfile <- SeqOpen(seqfile)
    snp.chr <- read.gdsn(index.gdsn(gfile, "snp.chromosome"))
    closefn.gds(gfile)
    rm(gfile)
    if(!"X" %in% snp.chr){
        message("\nDo not have chrX variants, skip sex check.\n")
        remove.samples <- remove.samples
    } else{
        seqfile <- SexCheck(seqfile)
        res.sexcheck <- QCresult(seqfile)$SexCheck
        table(res.sexcheck$sex, res.sexcheck$pred.sex)
        prob.sex <- res.sexcheck[res.sexcheck$sex != res.sexcheck$pred.sex & res.sexcheck$pred.sex != "0", ]
        remove.sex <- prob.sex[,1]
        remove.samples <- unique(c(remove.samples, remove.sex))
        if(results){
            write.table(res.sexcheck, file=paste0(dirname(output), "/result.sexcheck.txt"), quote=FALSE, sep="\t", row.names=FALSE)
        }
        if(plotting){
            p <- plotQC(seqfile, "SexCheck")
            ggsave(filename = paste0(dirname(output), "/plot.sexcheck.pdf"), p)
        }
    }
    
    ####################
    ## Inbreeding check
    
    seqfile <- Inbreeding(seqfile, remove.samples=remove.samples)
    res.inb <- QCresult(seqfile)$Inbreeding
    prob.inb <- res.inb[res.inb$outlier.5sd == "Yes" & !is.na(res.inb$outlier.5sd), ]
    remove.inb <- prob.inb[,1]
    remove.samples <- unique(c(remove.samples, remove.inb))
    if(results){
    write.table(res.inb, file=paste0(dirname(output), "/result.inbreeding.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    }   
    ## plot Inbreeding.
    if(plotting){
        p <- plotQC(seqfile, "Inbreeding")
        ggsave(filename = paste0(dirname(output), "/plot.inbreeding.pdf"), p)
    }
    
    ########
    ## IBD
    
    seqfile <- IBDCheck(seqfile, remove.samples=remove.samples)
    res.ibd <- QCresult(seqfile)$IBD
    prob.ibd <- IBDRemove(seqfile)
    ibd.pairs <- prob.ibd$ibd.pairs
    if(nrow(ibd.pairs) > 0){
        ibd.pairs <- paste(ibd.pairs$id1, ibd.pairs$id2, sep=":")
    }else{
        ibd.pairs = NULL
    }
    remove.ibd <- prob.ibd$ibd.remove
    remove.samples <- unique(c(remove.samples, remove.ibd))
    if(results){
        write.table(res.ibd, file=paste0(dirname(output), "/result.ibd.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    }
    ## plot IBD.
    if(plotting){
        p <- plotQC(seqfile, "IBD")
        ggsave(filename = paste0(dirname(output), "/plot.ibd.pdf"), p)
    }
    
    ########
    ## PCA

    seqfile <- PCACheck(seqfile, remove.samples=remove.samples)
    res.pca <- QCresult(seqfile)$PCA
    res.study <- res.pca[res.pca$sample %in% studyid, ]
    remove.pca <- res.study[res.study$pop != res.study$pred.pop, 1]

    if(results){
    write.table(res.pca, file=paste0(dirname(output), "/result.pca.txt"), quote=FALSE, sep="\t", row.names=FALSE)   
    }
    ## plot PCA.
    if(plotting){
        p <- plotQC(seqfile, "PCA")
        ggsave(filename = paste0(dirname(output), "/plot.pca.pdf"), p)
    }
    
    ################################
    ## summary of removed samples.
    problem.list <- problemList(seqfile)
    prob.list <- problem.list$prob.list
    rm.list <- problem.list$remove.list

    if(nrow(prob.list) != 0){
        a <- QCresult(seqfile)
        a$problem.list <- prob.list
        a$remove.list <- rm.list
        write.table(prob.list, file=paste0(dirname(output), "/result.problemSamples.txt"), quote=FALSE, sep="\t", row.names=FALSE)
        write.table(rm.list, file=paste0(dirname(output), "/result.removeSamples.txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }

    if(QCreport){
        RenderReport(seqfile, output=out.report, interactive=interactive)
    }
    
    ###########################
    ## save the seqfile in SeqSQC
    ###########################
    return(seqfile)
}
