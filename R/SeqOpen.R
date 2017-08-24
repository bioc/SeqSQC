#' Open the gds file with SeqSQC functions. 
#'
#' Function to open the gds file inside the SeqSQCclass file.
#' @param seqfile SeqSQCclass input file, which has been merged with benchmark data.
#' @param readonly whether to open the gds file in read-only mode. If "FALSE", it is allowed to write data to the file. The default is TRUE.
#' @param allow.duplicate whether to allow to open a GDS file with read-only mode when it has been opened in the same R session. The default is FALSE.
#' @export
#' @examples
#' library(gdsfmt)
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' seqfile@gdsfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' dat <- SeqOpen(seqfile)
#' dat
#' closefn.gds(dat)

SeqOpen <- function(seqfile, readonly=TRUE, allow.duplicate=FALSE)
{
    ## check
    fn <- seqfile@gdsfile
    stopifnot(is.character(fn), length(seqfile)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)

    ## open the file
    dat <- openfn.gds(fn, readonly=readonly, allow.fork=TRUE, allow.duplicate=allow.duplicate)
    return(dat)
}
