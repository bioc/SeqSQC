
subsetGDS <- function(gds, output.gds=NULL, sample.idx=NULL, snp.idx=NULL){

    ## when both "sample.idx" and "snp.idx" is defined, we remove snps first, and then remove samples, and then remove corresponding homozygous ref sites. 

    filepath <- gds$filename
    closefn.gds(gds)
    rm(gds)
    ## cleanup.gds(filepath)
    
    ## write to a new gds file. 
    if (!is.null(output.gds)){
        file.copy(filepath, output.gds, overwrite=TRUE)
        gds <- openfn.gds(output.gds, readonly=FALSE)
    }else{
        gds <- openfn.gds(filepath, readonly=FALSE)
    }
    
    if(!("gds.class" %in% class(gds)))stop("Not gds format!")
    nodes <- ls.gdsn(gds)
    snp.nodes <- nodes[grep("snp", nodes)]

    ## GDS subset for SNPs.
    if(!is.null(snp.idx)){
        for(i in snp.nodes){
            if(i=="snp.annot"){
                sa <- index.gdsn(gds, "snp.annot")
                for(j in ls.gdsn(sa)){
                    assign.gdsn(index.gdsn(sa, j), seldim = snp.idx)
                }
            }else {
                assign.gdsn(index.gdsn(gds, i), seldim = snp.idx)
            }
        }
        assign.gdsn(index.gdsn(gds, "genotype"), seldim = list(snp.idx, NULL))
    }
    
    ## GDS subset for samples, (do not check homozygous reference SNPs after removing samples.)
    sample.id <- index.gdsn(gds, "sample.id")
    if("sample.annot" %in% ls.gdsn(gds$root)){
        sample.annot <- index.gdsn(gds, "sample.annot")
    }else{
        sample.annot <- NULL
    }

    if(!is.null(sample.idx)){
        assign.gdsn(sample.id, seldim = sample.idx)
        assign.gdsn(index.gdsn(gds, "genotype"), seldim = list(NULL, sample.idx))
        if(!is.null(sample.annot)){
            for(i in ls.gdsn(sample.annot)){
                assign.gdsn(index.gdsn(sample.annot, i), seldim = sample.idx)
            }
        }
        
        ## ## check if any SNPs does not have variants after removing samples, and remove those SNPs.
        ## ## snp.idx.nhr <- unlist(apply.gdsn(index.gdsn(gds, "genotype"), 1, function(x) sum(x) > 0)) ## not correct, 0/1/2 is dosage of ref alleles, should count 2 for homo-ref here. 
        ## nsamples <- length(read.gdsn(sample.id))
        ## ## snp.idx.nhr <- unlist(apply.gdsn(index.gdsn(gds, "genotype"), 1, function(x) sum(x) < nsamples*2))
        ## snp.idx.nhr <- unlist(apply.gdsn(index.gdsn(gds, "genotype"), 1, function(x) sum(x==2) < nsamples))
        ## ## take ~20 min for whole genome seq data.
        
        ## if(any(!snp.idx.nhr)){
        ##     for(i in snp.nodes){
        ##         if(i=="snp.annot"){
        ##             sa <- index.gdsn(gds, "snp.annot")
        ##             for(j in ls.gdsn(sa)){
        ##                 assign.gdsn(index.gdsn(sa, j), seldim = snp.idx.nhr)
        ##             }
        ##         }else {
        ##             assign.gdsn(index.gdsn(gds, i), seldim = snp.idx.nhr)
        ##         }
        ##     }
        ##     assign.gdsn(index.gdsn(gds, "genotype"), seldim = list(snp.idx.nhr, NULL))
        ## }
    }

    closefn.gds(gds)
    if (!is.null(output.gds)){
        cleanup.gds(output.gds)
    }
    else{
        cleanup.gds(filepath)
    }
}
