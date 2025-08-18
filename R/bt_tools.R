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
#' Handles the following command:
#'  bedtools genomecov -ibam input_file -bg (-strand strand) > output_file.
#' Creates or returns coverage file.
#'
#'
#' @param input_bed Path to input BED file.
#' @param genome_file Path to genome file.
#' @param output_file [optional] path to output file.
#' The extension in the given name will be stripped.
#' If stranded = TRUE, then two output files with "_pos" and "_neg" appended
#' will be created.
#' @param stranded logical. Should output be split into strand-specific two files?
#' Uses -strand argument of genomecov.
#' Default is FALSE.
#'
#'
#' @return If output file given, its path or a vector of paths (if stranded) is returned.
#' Otherwise, a data table of the result captured from the command, or
#' if stranded, a list of two data tables: (+,-).
#' @export
#'
#' @examples
#' \dontrun{
#'   set_bedtools_path("/usr/bin/bedtools")
#'   bt_genomecov("bamfile.bam", "output.bedgraph", stranded = TRUE)
#' }
bt_genomecov <- function(input_bam, output_file = NULL, stranded = FALSE) {
  # Check input file type
  stopifnot(endsWith(input_bam, ".bam"))

  params <- as.list(environment())

  # ================ Config and Logs  ================
  write_config("make_windows", params)
  revert_sink <- start_log("bt_genomecov", params) # starts the log and returns a function to revert the sink on exit
  if (!is.null(revert_sink)) on.exit(revert_sink(), add = T)
  # ================================

  # 0: not stranded, 1: regular, -1: reverse stranded

  # stranded = FALSE: not stranded
  # stranded = TRUE: strand-specific output (split into + and -) -> two files

  bedtools <- get_bedtools_path()
  base_args <- c("genomecov", "-ibam", input_bam, "-bg")

  # Set strand mapping if given
  strand_map <- c(pos = "+", neg = "-")

  if (!is.null(output_file)) {
    output_file <- sub("\\.bedgraph$", "", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

    # if stranded -> split the result into strands
    if (stranded) {

      for (s in names(strand_map)) {
        args <- c(base_args, "-strand", strand_map[[s]])
        file <- paste0(output_file, "_", s, ".bedgraph")
        processx::run(bedtools, args, stdout = file, echo_cmd = T)
      }
    }

    # if not stranded do one call
    else {
      processx::run(bedtools, base_args, stdout = paste0(output_file, ".bedgraph"), echo_cmd = T)
    }

    return(invisible(output_file))
  }

  # No output file -> return data in R
  if (stranded) {
    result <- lapply(strand_map, function(s) {
      args <- c(base_args, "-strand", s)
      res <- processx::run(bedtools, args)
      data.table::fread(res$stdout)
    })
    names(result) <- names(strand_map)
    return(result)
  }

  res <- processx::run(bedtools, base_args)
  return(data.table::fread(res$stdout))
}

# HELPERS =========
