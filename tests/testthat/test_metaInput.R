context("sample meta info")

test_that("population input is correctks", {
    sample.annot <- system.file("extdata", "sampleAnnotation.txt", package="SeqSQC")
    meta <- read.table(sample.annot, header=T)
    expect_true(all(meta[,2] %in% c("AFR", "EUR", "EAS", "SAS", "ASN")))
})
