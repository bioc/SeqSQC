
mergeGDS <- function(gds1, gds2, output, missing.fill=TRUE){
    if(!("gds.class" %in% class(gds1)))stop("Not gds.class!")
    if(!("gds.class" %in% class(gds2)))stop("Not gds.class!")
    ## merge samples
        
    gds.out <- createfn.gds(output, allow.duplicate = TRUE)

    message("merging sample information...")
    ids <- c(read.gdsn(index.gdsn(gds1, "sample.id")),
             read.gdsn(index.gdsn(gds2, "sample.id")))
    anno <- rbind(read.gdsn(index.gdsn(gds1, "sample.annot")),
                  read.gdsn(index.gdsn(gds2, "sample.annot")))
    add.gdsn(gds.out, "sample.id", ids)
    add.gdsn(gds.out, "sample.annot", anno, check = FALSE)
    message("merging sample information... done!")
    
    ## merge snps
    message("merging snp information...")

    snpn <- c("snp.id", "snp.rs.id", "snp.chromosome", "snp.position",
              "snp.allele") 
    snpnodes1 <- lapply(
        snpn,
        function(x){
            read.gdsn(index.gdsn(gds1, x))
        }
    )
    snpnodes2 <- lapply(
        snpn,
        function(x){
            read.gdsn(index.gdsn(gds2, x))
        }
    )
    names(snpnodes1) <- names(snpnodes2) <- snpn

    snp1 <- paste(
        snpnodes1$snp.chromosome,
        snpnodes1$snp.position,
        snpnodes1$snp.allele,
        sep="_")
    snp2 <- paste(
        snpnodes2$snp.chromosome,
        snpnodes2$snp.position,
        snpnodes2$snp.allele,
        sep="_")
    
    if(missing.fill == TRUE){
        ## fill non-intersect genotypes as homozygous reference.
        d21 <- !(snp2 %in% snp1)
        d12 <- !(snp1 %in% snp2)
        snps1 <- c(snp1, snp2[d21])
        snps2 <- c(snp2, snp1[d12])
        
        add.gdsn(
            gds.out,
            "snp.id",
            c(snpnodes1$snp.id, max(snpnodes1$snp.id) + snpnodes2$snp.id[d21])
        )
        add.gdsn(gds.out,
                 "snp.rs.id",
                 c(snpnodes1$snp.rs.id, snpnodes2$snp.rs.id[d21])
                 )
        add.gdsn(
            gds.out,
            "snp.chromosome",
            c(snpnodes1$snp.chromosome, snpnodes2$snp.chromosome[d21])
        )
        add.gdsn(
            gds.out,
            "snp.position",
            c(snpnodes1$snp.position, snpnodes2$snp.position[d21])
        )
        add.gdsn(
            gds.out,
            "snp.allele",
            c(snpnodes1$snp.allele, snpnodes2$snp.allele[d21])
        )
        message("merging snp information... done!")
        
        message("merging genotypes...")
        add.gdsn(gds.out, "genotype", storage = "bit2",
                 valdim=c(length(snps1), 0))

        apply.gdsn(index.gdsn(gds1, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"),
                        c(x, rep(2, sum(d21)))) ## fill as homo ref.
        })
        apply.gdsn(index.gdsn(gds2, "genotype"), 2, function(x){
            append.gdsn(index.gdsn(gds.out, "genotype"),
                        c(x, rep(2, sum(d12)))[match(snps1, snps2)])
        })
        message("mergeing genotypes... done!")
        
        if("snp.annot" %in% ls.gdsn(gds1) & "snp.annot" %in% ls.gdsn(gds2)){
            message("merging snp annotation file for NGS data ...")
            
            snpanno <- c("snp.annot/qual", "snp.annot/filter")
            snpanno1 <- lapply(
                snpanno,
                function(x){
                    read.gdsn(index.gdsn(gds1, x))
                }
            )
            snpanno2 <- lapply(
                snpanno,
                function(x){
                    read.gdsn(index.gdsn(gds2, x))
                }
            )
            names(snpanno1) <- names(snpanno2) <- c("snp.qual", "snp.filter")
            
            af <- addfolder.gdsn(gds.out, "snp.annot")
            add.gdsn(af, "qual", c(snpanno1$snp.qual, snpanno2$snp.qual[d21]))
            add.gdsn(af, "filter", c(snpanno1$snp.filter, snpanno2$snp.filter[d21]))
            message("merging snp annotation file for NGS data ... done!")
        }
        
    }else{
        ## only keep the intersection snps.
        ind <- snp1 %in% snp2
        add.gdsn(gds.out, "snp.id", snpnodes1$snp.id[ind])
        add.gdsn(gds.out, "snp.rs.id", snpnodes1$snp.rs.id[ind])
        add.gdsn(gds.out, "snp.chromosome", snpnodes1$snp.chromosome[ind])
        add.gdsn(gds.out, "snp.position", snpnodes1$snp.position[ind])
        add.gdsn(gds.out, "snp.allele", snpnodes1$snp.allele[ind])
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

