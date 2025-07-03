.binGenome_env <- new.env(parent = emptyenv())
.binGenome_env$binGenome_path <- "binGenome.sh" # default, used if in the same dir

#' @rdname bin_genome
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
#' @param input_dir Path to the folder containing bedgraph files.
#' @param strand Filters files according to their strandedness. Currently a grep operation on "pos"/"neg".
#' @param annotation Annotation file (.bed)
#' @param output_dir Output folder for binGenome.sh
#' @param bins/quantiles bins or quantiles for binGenome
#' @param grep_pattern [optional] Filter for a grep pattern in files.
#' @param fixedBinSizeDownstream,fixedBinSizeUpstream [optional] Specifies the parameters that will be
#' used for binning the downstream or upstream regions, respectively. Format: "binsize:binnumber".
#' @param bedgraphNames,annotationNames [optional]
#' @param cores [optional] Number of threads to use. Defualt is max(1, parallel::detectCores() - 2)
#' @param normalize [optional]
#'
#'
#' @export
bin_genome <- function(input_dir, strand, annotation, output_dir,
                       bins = NULL, quantiles = NULL, grep_pattern = NULL,
                       fixedBinSizeDownstream = NULL, fixedBinSizeUpstream = NULL,
                       bedgraphNames = NULL, annotationNames = NULL,
                       cores = NULL, normalize = FALSE) {
  print(paste("Using programme:", .binGenome_env$binGenome_path))

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
  console_output <- processx::run(
    command = .binGenome_env$binGenome_path,
    args = args,
    echo = FALSE
  )
  too_small <- strsplit(console_output$stdout, "\n")[[1]] %>%
    base::grepl("[WARN].*(too small)", .) %>%
    sum
  print(paste("Too small regions found:", too_small))

}


