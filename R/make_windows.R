utils::globalVariables(c("gtf", "bed", "name", "path_to_out", "start", "end", "strand", "seqnames", "score", "gene_id", "input"))

#' Create genomic regions with upstream-downstream sizes and optional filters
#'
#' @param input_file GTF or BED file
#' @param upstream Upstream size
#' @param downstream Downstream size
#' @param path_to_output (Optional) File path to write output. Creates dirs if they do not exist
#' @param feature (Optional) Filter by given feature such as gene, transcript...
#' @param biotype (Optional) Filter by given biotype, if exists in the input GTF, such as lncRNA, pseudogene...
#' @param position (Optional) Specify 'start' or 'end' to build the window around.
#' If not specified, the start and end coordinates are extended with upstream and downstream sized in corresponding way.
#'
#' @returns A data.frame of annotated .bed file with specified windows and filter criteria.
#' The columns 'seqnames' and 'strand' are factors.
#' @export
#'
#' @examples
#' example_file <- system.file("extdata", "data", "example.gtf", package = "makeregions")
#' result <- make_windows(input_file=example_file,
#'                        upstream = 1000, downstream = 2000,
#'                        feature = "gene", biotype = "lncRNA",
#'                        position = "start")$result
#' print(head(result))
make_windows <- function(input_file, upstream, downstream,
                         path_to_output = NULL,
                         feature = NULL, biotype = NULL, position = NULL) {
  score <- NA
  # Argument checks
  if (missing(input_file) || missing(upstream) || missing(downstream)) {
    stop("Missing required arguments: --upstream, --downstream")
  }
  if(endsWith(input_file, ".gtf")) {
    input_file <- rtracklayer::import(input_file) %>%
      as.data.frame()
  } else if (endsWith(input_file, ".bed")) {
    input_file <- rtracklayer::import(input_file) %>%
      as.data.frame(strings) %>%
      dplyr::rename(gene_id = name)
  } else {
    stop("Provided input file is of wrong format. Should: .gtf or .bed")
  }

  # optional argument check
  if(!is.null(position) && !(position %in% c("start", "end"))) {
    stop("Position argument should be either 'start' or 'end' !")
  }
  # TODO: Check path_to_out validity

  # Create the directory path in out if they don't exist
  if (!is.null(path_to_output)) {
    dir_path <- dirname(path_to_out)
    if (dir_path != "." && !dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  # Save parameters
  params <- as.list(environment())

  # Create the .bed file
  annotation <- input_file %>%

    # Filter feature (gene, variaton, ...) and biotype (protein_coding, lncRNA, ...)
    # columns are quoted
    dplyr::filter(using_if_given(., "gene_biotype", biotype)) %>%
    dplyr::filter(using_if_given(., "type", feature)) %>%

    # if position argument is given, set both coordinates to the given position
    dplyr::mutate(start = if (!is.null(position) && position == "end") end else start,
                  end   = if (!is.null(position) && position == "start") start else end)%>%

    # Calculate upstream-downstream windows
    dplyr::mutate(start = ifelse(strand == "+", start-1-upstream, start-1-downstream),
           end = ifelse(strand == "+", end+downstream, end+upstream)) %>%

    # if score not found, fill with dot in .bed annotation
    dplyr::mutate(score = ifelse(is.na(score), ".", score)) %>%

    # remove regions where coordinates go over the edge
    dplyr::filter({
      negative_detected <- start < 0 | end < 0
      if (any(negative_detected, na.rm = TRUE)) {
        message(paste("Negative coordinates filtered out!", sum(negative_detected), "rows eliminated"))
      }
      !negative_detected
    }) %>%
    dplyr::select(seqnames, start, end, gene_id, score, strand)

  if(!is.null(path_to_output)) {
    annotation %>%
      readr::write_tsv(file.path(path_to_out), col_names = F)
  }


  return(structure(list(result = as.data.frame(annotation),
                        params = params),
                   class = "GenomicRegion"))
}



# HELPERS

# handle filter operations
using_if_given <- function(df, column_name, value) {
  # df:           the data frame
  # column_name:  the name of the column as a string (e.g., "biotype")
  # value:        the value to filter on (e.g., "protein_coding")

  # Ignore 1: If value not provided
  if (is.null(value) || is.na(value)) {
    cat("Argument not given:", deparse(substitute(value)), # backtracking for variable name
        "- continuing without filter.\n")
    return(TRUE)
  }

  # Ignore 2: If col does not exist
  if (!(column_name %in% colnames(df))) {
    cat("The input data does not have column:", column_name,
        "- continuing without filter.\n")
    return(TRUE)
  }

  # Success: Perform filter
  return(df[[column_name]] == value)
}


