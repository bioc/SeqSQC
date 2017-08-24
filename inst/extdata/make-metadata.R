library(AnnotationHubData)

meta <- data.frame(
    Title = "DNA-Sequencing dataset from the 1000 Genomes Project",
    Description = paste0("DNA-seq data from the 1000 Genomes Project ",
                         "containing 22 AFR, 22 EAS, 21 EUR and 22 SAS samples. ",
                         "there are eight known related pairs including ",
                         "four parent-offspring pairs, two full-sibling pairs ",
                         "and two half-sibling or avuncular pairs in this dataset, ",
                         "which is saved as a Genomic Data Structure (GDS) file."),
    BiocVersion = "3.4",
    Genome = "GRCh37",
    SourceType = "VCF",
    SourceUrl = paste0("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502",
                       "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/related_samples_vcf/", collapse=","),
    SourceVersion = "May 02 2013",
    Species = "Homo sapiens",
    TaxonomyId = 9606,
    Coordinate_1_based = TRUE,
    DataProvider = "1000 Genomes Project",
    Maintainer = "Qian Liu <qliu7@buffalo.edu>",
    RDataClass = "gds.class",
    DispatchClass = "GDSResource",
    ResourceName = "benchmark_1000genomes.gds"
)

## Not run:

## Write the data out and put in the inst/extdata directory.
write.csv(meta, file="seqsqc/inst/extdata/metadata.csv", row.names=FALSE)

## Test the validity of metadata.csv with readMetadataCsv():
readMetadataFromCsv("seqsqc")
## End(Not run)
