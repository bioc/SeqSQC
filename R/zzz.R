.onLoad <- function(libname = find.package("SeqSQC"), pkgname = "SeqSQC"){

    if(getRversion() >= "2.15.1"){
        utils::globalVariables(c("PC1", "PC2", "id1", "id2",
                                 "inbreeding", "k0", "k1", "label",
                                 "missingRate", "pop", "pred.label",
                                 "pred.pop", "pred.sex",
                                 "problematic", "relation", "sex",
                                 "sexinb", "type"))
    }
    invisible()
}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.20")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}

