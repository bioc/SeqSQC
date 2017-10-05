#' Data preprocessing for VCF or plink input from NGS or GWAS data.
#'
#' Function to read VCF or plink files, merge with benchmark data, and output as SeqSQCclass file.
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
#' @param ... Arguments to be passed to other methods.
#' @export
#' @return a SeqSQCclass object with the filepath to the gds file which stores the genotype, the summary of samples and variants, and the QCresults including the sample annotation information.  
#' @details
#' For \code{vfile} with more than one file names, \code{LoadVfile} will merge all dataset together if they all contain the same samples. It is useful to combine genetic/genomic data together if VCF data is divided by chromosomes. \cr
#' \code{sample.annot} file contains 3 columns with column names. col 1 is \code{sample} with sample ids; col 2 is \code{population} with values of "AFR/EUR/ASN/EAS/SAS"; col 3 is \code{gender} with values of "male/female".
#' @import gdsfmt
#' @import SNPRelate
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import methods
#' @import ExperimentHub
#' @importFrom utils read.table
#' @examples
#' infile <- system.file("extdata", "example_sub.vcf", package="SeqSQC")
#' sample.annot <- system.file("extdata", "sampleAnnotation.txt", package="SeqSQC")
#' cr <- system.file("extdata", "CCDS.Hs37.3.reduced_chr1.bed", package="SeqSQC")
#' outfile <- "testWrapUp"
#' seqfile <- LoadVfile(vfile = infile, output = outfile, capture.region = cr, sample.annot = sample.annot)
#' ## save(seqfile, file="seqfile.RData")
#' @author Qian Liu \email{qliu7@@buffalo.edu}

LoadVfile <- function(vfile, output, capture.region=NULL, sample.annot=NULL, LDprune=TRUE, vfile.restrict=FALSE, slide.max.bp=5e+05, ld.threshold=0.3, format.data="NGS", format.file="vcf", ...){

    tmpdir <- tempdir()
    study.gds <- tempfile(tmpdir=tmpdir)
    ## 1. read in vcf/plink file from study cohort, convert to gds format. 
    if(format.file == "vcf"){
        message("Load vcf file ...")
        ## study.gds <- "study.gds"
        snpgdsVCF2GDS(vfile, study.gds, method="biallelic.only", snpfirstdim=TRUE, ...)
    }else if(format.file == "plink"){
        message("Load plink file ...")
        ## study.gds <- "study.gds"
        vfile.base <- sub("ped|map|bed|bim|fam", "", vfile)
        if(length(grep("ped|map", vfile)) != 0){
            vfile.ped <- paste0(vfile.base, "ped")
            vfile.map <- paste0(vfile.base, "map")
            snpgdsPED2GDS(vfile.ped, vfile.map, study.gds, family=TRUE, snpfirstdim=TRUE, ...)
        }else{
            vfile.bed <- paste0(vfile.base, "bed")
            vfile.fam <- paste0(vfile.base, "fam")
            vfile.bim <- paste0(vfile.base, "bim")
            snpgdsBED2GDS(vfile.bed, vfile.fam, vfile.bim, study.gds, family=TRUE, snpfirstdim=TRUE, ...)
        }
    }else{
        stop("wrong file type")
0    }
    ## vtmp <- tempfile(tmpdir=tmpdir)
    ## file.copy(study.gds, vtmp)
    studycohort <- openfn.gds(study.gds, readonly=FALSE) 
    samples <- read.gdsn(index.gdsn(studycohort, "sample.id"))
    
    ## 2. read in the annotation for the study samples.
    if(!"sample.annot" %in% ls.gdsn(studycohort)){
        message("Load study cohort annotation file ...")
        study.annot <- read.table(file=sample.annot, as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)
        names(study.annot) <- c("sample", "population", "gender")
        study.annot$relation <- NA
        study.annot$group <- "study"
        study.annot <- study.annot[match(samples, study.annot$sample), ]
        add.gdsn(studycohort, "sample.annot", val=study.annot, check=FALSE)
    }else{
        message("Rewrite study cohort annotation file for PLINK data ...")
        study.annot <- read.table(file=sample.annot, as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)
        sampleanno <- read.gdsn(index.gdsn(studycohort, "sample.annot"))

        relations <- rep("", length(samples))
        ind1 <- sampleanno[, "father"] != 0 & sampleanno[, "mother"] != 0
        relations[ind1] <- paste0("father:", sampleanno[ind1, "father"], ";mother:", sampleanno[ind1, "mother"])
        ind2 <- sampleanno[, "father"] != 0 & sampleanno[, "mother"] == 0
        relations[ind2] <- paste0("father:", sampleanno[ind2, "father"])
        ind3 <- sampleanno[, "father"] == 0 & sampleanno[, "mother"] != 0
        relations[ind3] <- paste0("mother:", sampleanno[ind3, "mother"])

        ind <- match(samples, study.annot$sample)
        annot <- data.frame(sample=samples, population=study.annot[ind, "population"], gender=study.annot[ind, "gender"], relation=relations, group="study", stringsAsFactors=FALSE)
        add.gdsn(studycohort, "sample.annot", val=annot, check=FALSE, replace=TRUE)
    }        

    ## 3. read in benchmark gds file for 87 samples from 1000 genomes project.
    message("Load 1kg data to temp directory...")
    ## gds.1kg <- system.file("extdata", "benchmark_1000genomes.gds", package="SeqSQC")
    dfile <- ExperimentHub()[["EH550"]]
    gds.1kg <- dfile$filename
    closefn.gds(dfile)
    ## gds.1kg <- fileName(ExperimentHub())[["EH550"]]
    tmp.1kg <- tempfile(tmpdir=tmpdir)
    file.copy(gds.1kg, tmp.1kg)
    genofile <- openfn.gds(tmp.1kg, readonly = FALSE)

    ## 4. capture benchmark data (and/or) study data with Capture region
    if(!is.null(capture.region)){
        cp <- read.table(capture.region, sep="\t", stringsAsFactors=FALSE)
        cp <- GRanges(cp[,1], IRanges(as.numeric(cp[,2])+1, as.numeric(cp[,3])))
        
        chr.1kg <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
        pos.1kg <- read.gdsn(index.gdsn(genofile, "snp.position"))
        gr.1kg <- GRanges(chr.1kg, IRanges(pos.1kg, pos.1kg))
        ov <- findOverlaps(gr.1kg, cp)
        message("Subset 1kg data to capture region...")
        subsetGDS(genofile, snp.idx=unique(queryHits(ov)))
        genofile <- openfn.gds(tmp.1kg, readonly=TRUE)
    
        if(vfile.restrict==FALSE){
            chr.sc <- read.gdsn(index.gdsn(studycohort, "snp.chromosome"))
            pos.sc <- read.gdsn(index.gdsn(studycohort, "snp.position"))
            gr.sc <- GRanges(chr.sc, IRanges(pos.sc, pos.sc))
            message("restrict Vfile by capture region...")
            ov <- findOverlaps(gr.sc, cp)
            subsetGDS(studycohort, snp.idx = unique(queryHits(ov)))
            cleanup.gds(study.gds)
            studycohort <- openfn.gds(study.gds, readonly=TRUE)
        }
    }
    
    ## 5. merge study cohort with benchmark data
    message("Merging gds files of 1kg data and study cohort ...")
    output.merge <- paste0(output, ".gds")

    if(format.data == "NGS"){
        merge.out <- mergeGDS(gds1=genofile, gds2=studycohort, output=output.merge, missing.fill=TRUE)
        message("keep only snps within benchmark...")
        snp.bm <- paste(read.gdsn(index.gdsn(genofile, "snp.chromosome")), read.gdsn(index.gdsn(genofile, "snp.position")), read.gdsn(index.gdsn(genofile, "snp.allele")), sep="_")
        snp.merge <- paste(read.gdsn(index.gdsn(merge.out, "snp.chromosome")), read.gdsn(index.gdsn(merge.out, "snp.position")), read.gdsn(index.gdsn(merge.out, "snp.allele")), sep="_")
        bmtf <- snp.merge %in% snp.bm
        subsetGDS(gds=merge.out, snp.idx = bmtf)
        merge.out <- openfn.gds(output.merge, readonly=FALSE)
    }else{
        merge.out <- mergeGDS(gds1=genofile, gds2=studycohort, output=output.merge, missing.fill=FALSE)
    }
    class(merge.out) <- c("SNPGDSFileClass", "gds.class")
    
    ## LD pruning
    if(LDprune){
        message("LD pruning ...")
        ldprunein <- snpgdsLDpruning(merge.out, slide.max.bp=slide.max.bp, ld.threshold=ld.threshold)
        ldtf <- read.gdsn(index.gdsn(merge.out, "snp.id")) %in% unlist(ldprunein)
        if(! "snp.annot" %in% ls.gdsn(merge.out)){
            af <- addfolder.gdsn(merge.out, "snp.annot")
        }else{
            af <- index.gdsn(merge.out, "snp.annot")
        }
        add.gdsn(af, "LDprune", storage="logical", ldtf)
    }

    closefn.gds(genofile)
    closefn.gds(studycohort)
    unlink(tmp.1kg)
    unlink(study.gds)
    
    fn <- merge.out$filename
    samples <- read.gdsn(index.gdsn(merge.out, "sample.id"))
    snps <- read.gdsn(index.gdsn(merge.out, "snp.id"))
    sampleanno <- read.gdsn(index.gdsn(merge.out, "sample.annot")) 
    closefn.gds(merge.out)

    output <- SeqSQCclass(gdsfile = fn, QCresult = SimpleList(dimension = c(length(samples), length(snps)), sample.annot = sampleanno))
    return(output)
}
