#' Open the gds file in SeqSQC objects. 
#'
#' Function to open the gds file inside the SeqSQC object.
#' @param seqfile SeqSQC object, which has been merged with benchmark
#'     data.
#' @param readonly whether to open the gds file in read-only mode. If
#'     "FALSE", it is allowed to write data to the file. The default
#'     is TRUE.
#' @param allow.duplicate whether to allow to open a GDS file with
#'     read-only mode when it has been opened in the same R
#'     session. The default is FALSE.
#' @export
#' @return a gds file with the filepath in the input SeqSQC object.
#' @examples
#' library(gdsfmt)
#' load(system.file("extdata", "example.seqfile.Rdata", package="SeqSQC"))
#' gfile <- system.file("extdata", "example.gds", package="SeqSQC")
#' seqfile <- SeqSQC(gdsfile = gfile, QCresult = QCresult(seqfile))
#' dat <- SeqOpen(seqfile)
#' dat
#' closefn.gds(dat)
#' @author Qian Liu \email{qliu7@@buffalo.edu}

SeqOpen <- function(seqfile, readonly=TRUE, allow.duplicate=FALSE)
{
    ## check
    fn <- gdsfile(seqfile)
    stopifnot(is.character(fn), length(seqfile)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)

    ## open the file
    dat <- openfn.gds(fn, readonly=readonly, allow.fork=TRUE,
                      allow.duplicate=allow.duplicate)
    return(dat)
}
