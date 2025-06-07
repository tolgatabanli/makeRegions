
# makeRegions: Utility Tools to Create and Manipulate Genomic Regions

<!-- badges: start -->
<!-- badges: end -->

The package “makeRegions” currently helps with the creation of genomic
regions from GTF or BED files with options to filter specific features.

## Installation

You can install the development version of makeRegions from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("tolgatabanli/makeRegions")
```

## Example

``` r
library(makeRegions)

example_file <- system.file("extdata", "data", "example.gtf", package = "makeRegions")
result <- make_windows(input_file=example_file,
                       upstream = 1000, downstream = 2000,
                       feature = "gene", biotype = "lncRNA",
                       position = "start")$result
print(head(result))
#>   seqnames   start     end         gene_id score strand
#> 1     chr1 2580559 2583560 ENSG00000228037     .      +
#> 2     chr1 2579559 2582560 ENSG00000284616     .      -
#> 3     chr1 2580559 2583560 ENSG00000260972     .      +
#> 4     chr1 2579559 2582560 ENSG00000226374     .      -
#> 5     chr1 2580559 2583560 ENSG00000232596     .      +
```
