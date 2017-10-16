readgdsn <- function(gfile, node){
    if (!inherits(gfile, "gds.class")){
        return("object should inherit from 'gds.class'.")
    }
    stopifnot(class(node) == "character")
    stopifnot(node %in% ls.gdsn(gfile))

    gdsn <- read.gdsn(index.gdsn(gfile, node))
    gdsn
}
