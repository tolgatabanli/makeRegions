.binGenome_env <- new.env(parent = emptyenv())
.binGenome_env$binGenome_path <- "binGenome.sh" # default, used if in the same dir

#' Wrapper for binGenome module
#'
#' @rdname bin_genome
#' @name bin_genome
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

#' @rdname bin_genome
#' @description
#' Runs the binGenome.sh on the given file(s).
#' Prints the log of how many regions were too small for binning in each sample.
#' If stranded, pos and neg are merged and stripped from the name.
#'
#' @param input_dir Path to the folder containing bedgraph files.
#' Only the files with .bedgraph extension are considered.
#' @param strand Filters files according to their strandedness.
#' It applies grep and filters files having "pos"/"neg" anywhere in their names.
#' @param annotation Annotation file (.bed)
#' @param output_dir Output folder for binGenome.sh
#' @param bins/quantiles bins or quantiles for binGenome
#' @param grep_pattern [optional] Filter for a grep pattern in files.
#' @param fixedBinSizeDownstream,fixedBinSizeUpstream [optional] Specifies the parameters that will be
#' used for binning the downstream or upstream regions, respectively. Format: "binsize:binnumber".
#' @param bedgraphNames,annotationNames [optional]
#' @param cores [optional] Number of threads to use. Defualt: max(1, parallel::detectCores() - 2)
#' @param normalize [optional]
#'
#'
#' @export
bin_genome <- function(input_dir, strand, annotation, output_dir,
                       bins = NULL, quantiles = NULL, grep_pattern = NULL,
                       fixedBinSizeDownstream = NULL, fixedBinSizeUpstream = NULL,
                       bedgraphNames = NULL, annotationNames = NULL,
                       cores = NULL, normalize = FALSE) {
  message(paste("Using programme:", .binGenome_env$binGenome_path))


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

  print("Too small regions found:")
  print(too_small_warns_by_sample(console_output$stdout, strand))
}

### HELPERS ===============

too_small_warns_by_sample <- function(console_output, strand) {
  sample <- NULL
  log_by_sample <- list() # list will be named with samples/files
  log_lines <- strsplit(console_output, split = "\n")[[1]]

  for (line in log_lines) {
    # Start of a new file
    if (grepl("^\\[INFO\\] started with.*\\.bedgraph'", line)) {
      # Extract filename from path
      path_match <- regmatches(line, regexpr("'[^']+\\.bedgraph'", line))
      sample <- gsub(".*\\/|\\.bedgraph'", "", path_match)  # strip path and extension
      log_by_sample[[sample]] <- character()
    } else if (!is.null(sample)) { # Group the logs under current filename
      log_by_sample[[sample]] <- c(log_by_sample[[sample]], line)
    }
  }

  # log_by_sample is a list of vectors of lines
  too_smalls_by_sample <- sapply(log_by_sample,
                                 function(lines) {
                                   sum(base::grepl("[WARN].*(too small)", lines))
                                   })

  # if stranded merge neg and pos (one is actually always 0 but to be sure)
  if (strand != 0) {
    names(too_smalls_by_sample) <- gsub("(_pos|_neg)", "", names(too_smalls_by_sample))
    too_smalls_by_sample <- tapply(too_smalls_by_sample,
                                   names(too_smalls_by_sample),
                                   sum)
  }
  return(too_smalls_by_sample)
}
