#' The wrap-up function for sample QC of sequencing/GWAS data.
#' 
#' A wrap-up function for sample QC. It reads in the variant genotypes in vcf/PLINK format, merges study cohort with benchmark data, and performs sample QC for the merged dataset.
#' 
#' @param vfile vcf or PLINK input file (ped/map/bed/bim/fam with same basename). Vfile could be a vector of character strings, see details.
#' @param output a character string for name of merged data in SeqSQCclass. 
#' @param capture.region the BED file of sequencing capture regions. The default is NULL. For exome-sequencing data, the capture region file must be provided.
#' @param sample.annot sample annotation file with 3 columns including the sample id, sample population and sex info. The default is NULL.
#' @param LDprune whether to use LD-pruned snp set. The default is TRUE.
#' @param vfile.restrict whether the input vcf or plink file has already been restricted by capture region. The default is FALSE.
#' @param slide.max.bp the window size of SNPs when calculating linkage disequilibrium. The default is 5e+05. 
#' @param ld.threshold the r^2 threshold for LD-based SNP pruning if \code{LDprune = TRUE}. The default is 0.3.
#' @param format.data the data source. The default is \code{NGS} for sequencing data. 
#' @param format.file the data format. The default is \code{vcf}.
#' @param QCreport Whether to generate the sample QC report in html format.
#' @param out.report the file name for the sample QC report. The default is \code{report.html}.
#' @param interactive whether to generate interactive plots in the sample QC report if \code{QCreport = TRUE}.
#' @param ... Arguments to be passed to other methods.
#' @export
#' @return a SeqSQCclass object with the filepath to the gds file which stores the genotype, the summary of samples and variants, and the QCresults including the sample annotation information and all QC results.  

#' @details
#' For \code{vfile} with more than one file names, \code{sampleQC} will merge all dataset together if they all contain the same samples. It is useful to combine genetic/genomic data together if VCF data is divided by chromosomes. \cr
#' There are 3 columns in \code{sample.annot} file. col 1 is \code{sample} with sample ids, col 2 is \code{population} with values of "AFR/EUR/ASN/EAS/SAS", col 3 is \code{gender} with values of "male/female".
#' @importFrom utils write.table
#' @examples
#' \dontrun{
#' seqfile <- sampleQC(vfile = system.file("extdata", "example_sub.vcf", package="SeqSQC"),
#' output = "testWrapUp",
#' capture.region = system.file("extdata", "CCDS.Hs37.3.reduced.bed", package="SeqSQC"),
#' sample.annot = system.file("extdata", "sampleAnnotation.txt", package="SeqSQC"),
#' QCreport = TRUE, interactive = TRUE)
#' save(seqfile, "seqfile.RData")
#' }

sampleQC <- function(vfile, output, capture.region = NULL, sample.annot = NULL, LDprune = TRUE, vfile.restrict = FALSE, slide.max.bp = 5e+05, ld.threshold = 0.3, format.data = "NGS", format.file = "vcf", QCreport = TRUE, out.report="report.html", interactive = TRUE, ...){

    seqfile <- LoadVfile(vfile = vfile, output = output, capture.region = capture.region, sample.annot = sample.annot, ...)

    fn <- seqfile@gdsfile
    print(paste("gds file generated:", fn))

    sampleanno <- seqfile@QCresult$sample.annot
    samples <- sampleanno$sample
    studyid <- sampleanno[sampleanno[,5] == "study", 1]
    
    ##############
    ## sample missing rate
    
    seqfile <- MissingRate(seqfile)
    res.mr <- seqfile@QCresult$MissingRate
    res.study <- res.mr[res.mr$sample %in% studyid, ]
    prob.mr <- res.study[res.study$outlier == "Yes", ]
    remove.mr <- prob.mr[,1]
    remove.samples <- unique(remove.mr)
    write.table(res.mr, file=paste0(dirname(output), "/result.missingrate.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    p <- plotMissingRate(seqfile)
    ggsave(filename = paste0(dirname(output), "/plot.missingrate.pdf"), p)

    ##############
    ## Sex Check
    
    seqfile <- SexCheck(seqfile)
    res.sexcheck <- seqfile@QCresult$SexCheck
    table(res.sexcheck$sex, res.sexcheck$pred.sex)
    prob.sex <- res.sexcheck[res.sexcheck$sex != res.sexcheck$pred.sex & res.sexcheck$pred.sex != "0", ]
    remove.sex <- prob.sex[,1]
    remove.samples <- unique(c(remove.samples, remove.sex))
    write.table(res.sexcheck, file=paste0(dirname(output), "/result.sexcheck.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    
    p <- plotSexCheck(seqfile)
    ggsave(filename = paste0(dirname(output), "/plot.sexcheck.pdf"), p)

    ####################
    ## Inbreeding check
    
    seqfile <- Inbreeding(seqfile, remove.samples=remove.samples)
    res.inb <- seqfile@QCresult$Inbreeding
    prob.inb <- res.inb[res.inb$outlier.5sd == "Yes" & !is.na(res.inb$outlier.5sd), ]
    remove.inb <- prob.inb[,1]
    remove.samples <- unique(c(remove.samples, remove.inb))
    write.table(res.inb, file=paste0(dirname(output), "/result.inbreeding.txt"), quote=FALSE, sep="\t", row.names=FALSE)
        
    ## plot Inbreeding.
    p <- plotInbreeding(seqfile)
    ggsave(filename = paste0(dirname(output), "/plot.inbreeding.pdf"), p)

    
    ########
    ## IBD
    
    seqfile <- IBDCheck(seqfile, remove.samples=remove.samples)
    res.ibd <- seqfile@QCresult$IBD
    prob.ibd <- IBDRemove(seqfile)
    ibd.pairs <- prob.ibd$ibd.pairs
    if(nrow(ibd.pairs) > 0){
        ibd.pairs <- paste(ibd.pairs$id1, ibd.pairs$id2, sep=":")
    }else{
        ibd.pairs = NULL
    }
    remove.ibd <- prob.ibd$ibd.remove
    remove.samples <- unique(c(remove.samples, remove.ibd))
    write.table(res.ibd, file=paste0(dirname(output), "/result.ibd.txt"), quote=FALSE, sep="\t", row.names=FALSE)
        
    ## plot IBD.
    p <- plotIBD(seqfile)
    ggsave(filename = paste0(dirname(output), "/plot.ibd.pdf"), p)

    ########
    ## PCA

    seqfile <- PCACheck(seqfile, remove.samples=remove.samples)
    res.pca <- seqfile@QCresult$PCA
    res.study <- res.pca[res.pca$sample %in% studyid, ]
    remove.pca <- res.study[res.study$pop != res.study$pred.pop, 1]

    write.table(res.pca, file=paste0(dirname(output), "/result.pca.txt"), quote=FALSE, sep="\t", row.names=FALSE)   
    ## plot PCA.
    p <- plotPop(seqfile)
    ggsave(filename = paste0(dirname(output), "/plot.pca.pdf"), p)

    ################################
    ## summary of removed samples.
    prob.list <- problemList(seqfile)
    rm.list <- prob.list$remove
    ## prob.list <- data.frame(sample=c(remove.mr, remove.sex, remove.inb, ibd.pairs, remove.pca), remove.reason=c(rep("high missing rate", length(remove.mr)), rep("gender mismatch", length(remove.sex)), rep("inbreeding outlier", length(remove.inb)), rep("cryptic relationship", length(ibd.pairs)), rep("population outlier", length(remove.pca))))
    ## rm.list <- data.frame(sample = c(remove.mr, remove.sex, remove.inb, remove.ibd, remove.pca))

    if(nrow(prob.list) != 0){
        seqfile@QCresult$problem.list <- prob.list
        seqfile@QCresult$remove.list <- rm.list
        write.table(prob.list, file=paste0(dirname(output), "/result.problemSamples.txt"), quote=FALSE, sep="\t", row.names=FALSE)
        write.table(rm.list, file=paste0(dirname(output), "/result.removeSamples.txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }

    if(QCreport){
        RenderReport(seqfile, output=out.report, interactive=interactive)
    }
    
    ###########################
    ## save the seqfile in SeqSQCclass
    ###########################
    return(seqfile)
}
