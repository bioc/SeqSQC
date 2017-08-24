
mergeGDS <- function(gds1, gds2, output, missing.fill=TRUE){
    if(!("gds.class" %in% class(gds1)))stop("Not gds.class!")
    if(!("gds.class" %in% class(gds2)))stop("Not gds.class!")
    ## merge samples
        
    gds.out <- createfn.gds(output, allow.duplicate = TRUE)

    message("merging sample information...")
    ids <- c(read.gdsn(index.gdsn(gds1, "sample.id")), read.gdsn(index.gdsn(gds2, "sample.id")))
    anno <- rbind(read.gdsn(index.gdsn(gds1, "sample.annot")), read.gdsn(index.gdsn(gds2, "sample.annot")))
    add.gdsn(gds.out, "sample.id", ids)
    add.gdsn(gds.out, "sample.annot", anno)
    message("merging sample information... done!")
    
    ## merge snps
    message("merging snp information...")
    snp.id.1 <- read.gdsn(index.gdsn(gds1, "snp.id"))
    snp.rs.id.1 <- read.gdsn(index.gdsn(gds1, "snp.rs.id"))
    snp.chromosome.1 <- read.gdsn(index.gdsn(gds1, "snp.chromosome"))
    snp.position.1 <- read.gdsn(index.gdsn(gds1, "snp.position"))
    snp.allele.1 <- read.gdsn(index.gdsn(gds1, "snp.allele"))
    snp1 <- paste(snp.chromosome.1, snp.position.1, snp.allele.1, sep="_")
    ## snp11 <- paste(snp.chromosome.1, snp.position.1, sep="_")
    
    snp.id.2 <- read.gdsn(index.gdsn(gds2, "snp.id"))
    snp.rs.id.2 <- read.gdsn(index.gdsn(gds2, "snp.rs.id"))
    snp.chromosome.2 <- read.gdsn(index.gdsn(gds2, "snp.chromosome"))
    snp.position.2 <- read.gdsn(index.gdsn(gds2, "snp.position"))
    snp.allele.2 <- read.gdsn(index.gdsn(gds2, "snp.allele"))
    snp2 <- paste(snp.chromosome.2, snp.position.2, snp.allele.2, sep="_")
    ## snp22 <- paste(snp.chromosome.2, snp.position.2, sep="_")
    
    if(missing.fill == TRUE){
        ## fill non-intersect genotypes as homozygous reference.
        d21 <- !(snp2 %in% snp1)
        d12 <- !(snp1 %in% snp2)
        snps1 <- c(snp1, snp2[d21])
        snps2 <- c(snp2, snp1[d12])
        
        add.gdsn(gds.out, "snp.id", c(snp.id.1, max(snp.id.1) + snp.id.2[d21]))
        add.gdsn(gds.out, "snp.rs.id", c(snp.rs.id.1, snp.rs.id.2[d21]))
        add.gdsn(gds.out, "snp.chromosome", c(snp.chromosome.1, snp.chromosome.2[d21]))
        add.gdsn(gds.out, "snp.position", c(snp.position.1, snp.position.2[d21]))
        add.gdsn(gds.out, "snp.allele", c(snp.allele.1, snp.allele.2[d21]))
        message("merging snp information... done!")
        
        message("merging genotypes...")
        add.gdsn(gds.out, "genotype", storage = "bit2", valdim=c(length(snps1), 0))

        apply.gdsn(index.gdsn(gds1, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"), c(x, rep(2, sum(d21)))) ## fill as homo ref.
        })
        apply.gdsn(index.gdsn(gds2, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"), c(x, rep(2, sum(d12)))[match(snps1, snps2)])
        })
        message("mergeing genotypes... done!")
        
        if("snp.annot" %in% ls.gdsn(gds1) & "snp.annot" %in% ls.gdsn(gds2)){
            message("merging snp annotation file for NGS data ...")
            snp.qual.1 <- read.gdsn(index.gdsn(gds1, "snp.annot/qual"))
            snp.filter.1 <- read.gdsn(index.gdsn(gds1, "snp.annot/filter"))
            snp.qual.2 <- read.gdsn(index.gdsn(gds2, "snp.annot/qual"))
            snp.filter.2 <- read.gdsn(index.gdsn(gds2, "snp.annot/filter"))
            af <- addfolder.gdsn(gds.out, "snp.annot")
            add.gdsn(af, "qual", c(snp.qual.1, snp.qual.2[d21]))
            add.gdsn(af, "filter", c(snp.filter.1, snp.filter.2[d21]))
            message("merging snp annotation file for NGS data ... done!")
        }
        
    }else{
        ## only keep the intersection snps.
        ind <- snp1 %in% snp2
        add.gdsn(gds.out, "snp.id", snp.id.1[ind])
        add.gdsn(gds.out, "snp.rs.id", snp.rs.id.1[ind])
        add.gdsn(gds.out, "snp.chromosome", snp.chromosome.1[ind])
        add.gdsn(gds.out, "snp.position", snp.position.1[ind])
        add.gdsn(gds.out, "snp.allele", snp.allele.1[ind])
        message("merging snp information... done!")
        
        message("merging genotypes...")
        add.gdsn(gds.out, "genotype", storage = "bit2", valdim=c(sum(ind), 0))
        ind2 <- match(snp1, snp2)
        ind22 <- ind2[!is.na(ind2)]
        apply.gdsn(index.gdsn(gds1, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"), x[ind])
        })
        apply.gdsn(index.gdsn(gds2, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"), x[ind22])
        })
        message("mergeing genotypes... done!")
        ## message("filtering variants with MAF<0.05...")
        ## mr <- apply.gdsn(index.gdsn(gds.out, "genotype"), 1, function(x){  
        ##     1 - sum(x %in% c(0, 1, 2))/length(x)
        ## })
    }
    return(gds.out)
}

