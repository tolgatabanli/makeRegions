.bedtools_env <- new.env(parent = emptyenv())
.bedtools_env$bedtools_path <- "bedtools" # default, used if in PATH

#' bedtools: utilities for genomic regions
#'
#' @name bedtools
#' @title bedtools: utilities for genomic regions
#' @rdname bedtools
NULL

#' @rdname bedtools
#' @section set_bedtools_path:
#'
#' @param path Path to the bedtools executable.
#'
#' @return TRUE on successful execution, FALSE otherwise.
#' @export
set_bedtools_path <- function(path) {
  if (!file.exists(path)) stop("Invalid bedtools path: ", path)
  .bedtools_env$bedtools_path <- normalizePath(path)
  invisible(TRUE)
}

#' @rdname bedtools
#' @section get_bedtools_path:
#' @export
get_bedtools_path <- function() {
  path <- .bedtools_env$bedtools_path
  if (is.null(path)) stop("Bedtools path not set. Please call set_bedtools_path().")
  path
}

#' Get the genome coverage using bedtools
#'
#' @description
#' Handles following command: bedtools genomecov -ibam input_file -bg > output_file.
#' Used for its side effect, creating coverage file.
#'
#'
#' @param input_bed Path to input BED file.
#' @param genome_file Path to genome file.
#' @param output_file Optional path to output file.
#'
#' @return If output file path given, a string of its path is returned,
#' otherwise a data table of the result captured from the command is returned.
#' @export
#'
#' @examples
#' \dontrun{
#'   set_bedtools_path("/usr/bin/bedtools")
#'   bt_genomecov("bamfile.bam", "output.bedgraph")
#' }
bt_genomecov <- function(input_bam, output_file = NULL) {
  # Check input file type
  stopifnot(!endsWith(input_bam, "bam"))

  bedtools <- get_bedtools_path()
  args <- c("genomecov", "-ibam", input_bam, "-bg")

  if (!is.null(output_file)) {
    processx::run(bedtools, args, stdout = output_file)
    return(invisible(output_file))
  } else {
    res <- processx::run(bedtools, args)
    # stdout is a tab separated string, should parse with fread
    return(data.table::fread(res$stdout))
  }
}
