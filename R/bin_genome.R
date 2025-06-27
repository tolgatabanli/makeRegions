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
bin_genome <- function(annotation, bedgraph = NULL, bedgraphNeg = NULL, bedgraphPos = NULL,
                       outputDir,
                       bins = NULL, quantiles = NULL,
                       fixedBinSizeDownstream = NULL, fixedBinSizeUpstream = NULL,
                       bedgraphNames = NULL, annotationNames = NULL,
                       cores = NULL, normalize = FALSE) {
  # exclusive arguments
  if (!is.null(bedgraph) && (!is.null(bedgraphNeg) || !is.null(bedgraphPos))) {
    stop("Use either `bedgraph` OR both `bedgraphNeg` and `bedgraphPos`, not both.")
  }
  if (is.null(bedgraph) && (is.null(bedgraphNeg) || is.null(bedgraphPos))) {
    stop("You must provide either `bedgraph`, OR both `bedgraphNeg` and `bedgraphPos`.")
  }
  if (!is.null(bins) && !is.null(quantiles)) {
    stop("Use either `bins` OR `quantiles`, not both.")
  }
  if (is.null(bins) && is.null(quantiles)) {
    stop("You must provide either `bins` or `quantiles`.")
  }

  # save arguments in a named list
  opts <- list(
    "--annotation" = annotation,
    "--bedgraph" = bedgraph,
    "--bedgraphNeg" = bedgraphNeg,
    "--bedgraphPos" = bedgraphPos,
    "--outputDir" = outputDir,
    "--bins" = bins,
    "--quantiles" = quantiles,
    "--fixedBinSizeDownstream" = fixedBinSizeDownstream,
    "--fixedBinSizeUpstream" = fixedBinSizeUpstream,
    "--bedgraphNames" = bedgraphNames,
    "--annotationNames" = annotationNames,
    "--cores" = cores
    )

  # flatten -> c("--flag", "val", ...)
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
bin_genome.multiple(input_folder, strand = NULL, grep_pattern = NULL) {
  # strand = c("none", "pos", "neg", "group")
  files <- list.files(input_folder)

  # Take bedgraphs
  extension <- ".bedgraph"
  files <- files[grepl(paste0("\\.", extension, "$"), files)]

  if (!is.null(pattern)) {
    files <- grep(pattern, files, value = TRUE)
  }
}
