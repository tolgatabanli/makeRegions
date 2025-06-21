
# makeRegions: Utility Tools to Create and Manipulate Genomic Regions

<!-- badges: start -->
<!-- badges: end -->

The package “makeRegions” currently helps with the creation of genomic
regions from GTF or BED files with options to filter specific features.

## Installation

You can install the development version of makeRegions from
[GitHub](https://github.com/) with: (If cannot use any packages for installing from GitHub, see below)

```r
# install.packages("pak")
pak::pak("tolgatabanli/makeRegions")
```
or 
```r
# install.packages("devtools")
devtools::install_github("tolgatabanli/makeRegions")
```
If you need an alternative to use this package as a user- or directory-specific library (e.g. for a remote cluster) without devtools' install_github, pak or similar:

For Linux:
1. Clone the repository locally (where you manage the packages on your own).
2. In RStudio, make sure you are at the repository root.
3. Use ```devtools:build()```. This will generate a tar.gz file.
4. Scp the tar.gz to the remote server.
5. Use ```install.packages("PACKAGENAME.tar.gz", repo = NULL, lib="YOUR_LOCAL_LIB")```
6. Whenever you need to attach the package, use ```library("PACKAGENAME", lib.loc="path/to/YOUR_LOCAL_LIB")```.

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


