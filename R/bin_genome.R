.binGenome_env <- new.env(parent = emptyenv())
.binGenome_env$binGenome_path <- "binGenome.sh" # default, used if in the same dir

#' @rdname binGenome
#' @section set_binGenome_path:
#'
#' @param path Path to the binGenome.sh executable.
#'
#' @return TRUE on successful execution, FALSE otherwise.
#' @export
set_binGenome_path <- function(path) {
  if (!file.exists(path)) stop("Invalid binGenome.sh path: ", path)
  .binGenome_env$binGenome_path <- normalizePath(path)
  invisible(TRUE)
}


#' Wrapper for binGenome module
#'
#' @description
#' Runs the binGenome.sh on the given file(s).
#'
#'
#' @returns
#' @export
#'
#' @examples
bin_genome <- function(input_dir, strand, annotation, output_dir,
                       bins = NULL, quantiles = NULL, grep_pattern = NULL,
                       fixedBinSizeDownstream = NULL, fixedBinSizeUpstream = NULL,
                       bedgraphNames = NULL, annotationNames = NULL,
                       cores = NULL, normalize = FALSE) {
  # mutually exclusive parameters
  if (!is.null(bins) && !is.null(quantiles)) {
    stop("Use either `bins` OR `quantiles`, not both.")
  }
  if (is.null(bins) && is.null(quantiles)) {
    stop("You must provide either `bins` or `quantiles`.")
  }
  # 0: not stranded, 1: regular, -1: reverse stranded
  if (!(strand %in% c(-1, 0, 1))) {
    stop(paste("Unrecognized argument for strand:", strand,
               "\nShould be one of c(-1, 0, 1)"))
  }

  # Parse files
  files <- list.files(input_dir)

  # Take bedgraphs if other files exist
  extension <- "bedgraph"
  files <- grep(paste0("\\.", extension, "$"), files, value = T)

  if (!is.null(grep_pattern)) {
    files <- grep(grep_pattern, files, value = TRUE)
  }

  # If stranded, group files
  pos <- neg <- bedgraph <- bedgraphPos <- bedgraphNeg <- NULL
  if (strand != 0) {
    if (strand == -1) {
      pos <- "neg"
      neg <- "pos"
    } else if (strand == 1) {
      pos <- "pos"
      neg <- "neg"
    }
    bedgraphPos <- paste(file.path(input_dir, grep(pos, files, value = T)), collapse = ",")
    bedgraphNeg <- paste(file.path(input_dir, grep(neg, files, value = T)), collapse = ",")
  } else {
    bedgraph <- files
  }

  # if cores not given, automatically detect
  cores <- max(1, parallel::detectCores() - 2)

  # Create the out directory if it doesn't exist
  out_path <- dirname(output_dir)
  if (out_path != "." && !dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }


  # save arguments in a named list, NULL is okay
  opts <- list(
    "--annotation" = annotation,
    "--bedgraph" = bedgraph,
    "--bedgraphNeg" = bedgraphNeg,
    "--bedgraphPos" = bedgraphPos,
    "--outputDir" = output_dir,
    "--bins" = bins,
    "--quantiles" = quantiles,
    "--fixedBinSizeDownstream" = fixedBinSizeDownstream,
    "--fixedBinSizeUpstream" = fixedBinSizeUpstream,
    "--bedgraphNames" = bedgraphNames,
    "--annotationNames" = annotationNames,
    "--cores" = cores
    )

  # flatten -> c("--flag", "val", ...) and removes the NULL args
  args <- unlist(Map(
    function(flag, val) {
      if (!is.null(val)) c(flag, as.character(val)) else NULL
      }, names(opts), opts), use.names = FALSE)
  # add normalize if given
  if (normalize) {
    args <- c("--normalize")
  }

  # Run process
  processx::run(
    command = .binGenome_env$binGenome_path,
    args = args,
    echo = TRUE
  )
}



# TODO: take input folder, optionally filter by pos/neg, or specific grep
# TODO: take argument for reverse, unstranded, regular
# TODO: log info on how many regions too small

