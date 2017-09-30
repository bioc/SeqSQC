#' Render the rmarkdown file for generating a sample QC report.
#' 
#' Function to render the pre-compiled rmarkdown file to generate the sample QC report. 
#' @param input SeqSQCclass object with QC results.
#' @param output a character string to define the file name for the QC report.
#' @param interactive whether to generate interative plots in the report. The default is TRUE.
#' @importFrom rmarkdown render
#' @return Will incure the rendering of the rmarkdown file for generating the sample QC report. The report will return to the file denoted in \code{output} in the function. 
#' @export 
#' @author Qian Liu \email{qliu7@@buffalo.edu}

RenderReport <- function(input, output, interactive = TRUE){
    report <- system.file("extdata", "sampleQCReport.Rmd", package="SeqSQC")
    render(report,
           output_file = basename(output),
           output_dir=dirname(output),
           params=list(
               seqfile = input,
               interactive = interactive
           ))
}
