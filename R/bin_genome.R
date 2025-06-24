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
create_bins <- function(bedgraph, bedgraphPos, bedgraphNeg, annotation,
                        bins, quantiles, outputDir,
                        bedgraphNames = NULL, annotationNames = NULL,
                        cores = NULL, normalize) {
  # TODO: take input folder, optionally filter by pos/neg, or specific grep
  # TODO: take argument for reverse, unstranded, regular
  # TODO: log info on how many regions too small
  # TODO: use binGenome.sh

  # Define input and output directories
  inputdir <- "/path/to/input"

  # Construct full command
  res <- processx::run(
    command = "binGenome.sh",
    args = c(
      "--bedgraphNeg",
      paste(file.path(inputdir, "12AS_Ctl_1_pos.bedgraph"),
            file.path(inputdir, "12AS_Ctl_2_pos.bedgraph"), sep = ","),
      "--bedgraphPos",
      paste(file.path(inputdir, "12AS_Ctl_1_neg.bedgraph"),
            file.path(inputdir, "12AS_Ctl_2_neg.bedgraph"), sep = ","),
      "--annotation", "anno_100Up_200Down.bed",
      "--bins", "100",
      "--fixedBinSizeDownstream", "40:50",
      "--fixedBinSizeUpstream", "50:2",
      "--outputDir", outputDir
    ),
    echo = TRUE  # optional: shows output live in console
  )

}

