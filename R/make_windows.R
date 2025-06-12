utils::globalVariables(c("gtf", "bed", "name", "path_to_out", "start", "end", "strand", "seqnames", "score", "gene_id", "input", ":="))

#' Create genomic regions with upstream-downstream sizes and optional filters
#'
#' @param input_file GTF or BED file
#' @param upstream Upstream size
#' @param downstream Downstream size
#' @param path_to_output (Optional) File path to write output. Creates dirs if they do not exist
#' @param position (Optional) Specify 'start' or 'end' to build the window around.
#' If not specified, the start and end coordinates are extended with upstream and downstream sizes in corresponding way.
#' @param ... (Optional) Arguments to filter. For example, type (gene, transcript)
#' gene_biotype (lncRNA, pseudogene)...
#' If for any column a vector of length 2 or more is given, the regions will be generated
#' that satisfy ANY one of the elements.
#'
#' @returns A data.frame in .bed format with specified window and filter criteria.
#' The columns 'seqnames' and 'strand' are factors.
#' @export
#'
#' @examples
#' example_file <- system.file("extdata", "data", "example.gtf", package = "makeRegions")
#' result <- make_windows(input_file=example_file,
#'                        upstream = 1000, downstream = 2000,
#'                        position = "start",
#'                        type = "gene", gene_biotype = "lncRNA")$result
#' print(head(result))
make_windows <- function(input_file, upstream, downstream,
                         path_to_output = NULL, position = NULL,
                         ...) {
  score <- NA
  gtf_filters <- list(...)

  # Argument checks
  if (missing(input_file) || missing(upstream) || missing(downstream)) {
    stop("Missing required arguments: upstream, downstream")
  }
  if(endsWith(input_file, ".gtf")) {
    input_file <- rtracklayer::import(input_file) %>%
      as.data.frame()
  } else if (endsWith(input_file, ".bed")) {
    input_file <- rtracklayer::import(input_file) %>%
      as.data.frame() %>%
      dplyr::rename(gene_id = name)
    if(length(gtf_filters) == 0) {
      stop("GTF filters cannot be used with input file of format BED.")
    }
  } else {
    stop(paste0("Provided input file (", input_file, ") is of wrong format. Should: .gtf or .bed"))
  }

  # optional argument check
  if(!is.null(position) && !(position %in% c("start", "end"))) {
    stop("Position argument should be either 'start' or 'end' !")
  }

  # TODO: Check path_to_out validity

  # Create the directory path in out if they don't exist
  if (!is.null(path_to_output)) {
    dir_path <- dirname(path_to_output)
    if (dir_path != "." && !dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  # Save parameters
  params <- as.list(environment())

  #  ========================== Creating BED ==========================
  annotation <- data.table::as.data.table(input_file)

  # filter by gtf arguments (feature and 9th column of gtf)
  for (key in names(gtf_filters)) {
    val <- gtf_filters[[key]]

    if (is.null(val) || all(is.na(val))) {
      message("Argument not given for ", key, " — skipping filter.")
      next
    }

    if (!(key %in% names(annotation))) {
      message("Column ", key, " not found — skipping filter.")
      next
    }

    annotation <- annotation[get(key) %in% val]
  }

  # Adjust start and end based on position
  if (!is.null(position) && position == "end") {
    annotation[, start := end]
  } else if (!is.null(position) && position == "start") {
    annotation[, end := start]
  }

  # Apply upstream/downstream logic
  annotation[, `:=`(
    start = ifelse(strand == "+", start - 1 - upstream, start - 1 - downstream),
    end   = ifelse(strand == "+", end + downstream, end + upstream)
  )]

  # if score not found, fill with dot in .bed annotation
  annotation[, score := as.character(score)]
  annotation[is.na(score), score := "."]

  # Remove regions where coordinates go over the edge with a messaage
  negative_detected <- annotation[, start < 0 | end < 0]
  if (any(negative_detected, na.rm = TRUE)) {
    message(paste("Negative coordinates filtered out!", sum(negative_detected, na.rm = TRUE), "rows eliminated"))
  }
  annotation <- annotation[!negative_detected]

  annotation <- annotation[, .(seqnames, start, end, gene_id, score, strand)]


  # row number for tracking changes and logging
  message(paste("Input had", nrow(input_file),
                "rows, output has", nrow(annotation), "rows.",
                "Output has", data.table::uniqueN(annotation), "unique rows."))

  if(!is.null(path_to_output)) {
    annotation %>%
      readr::write_tsv(file.path(path_to_output), col_names = F)
    message(paste("File has been written:", path_to_output))
  }


  return(structure(list(result = as.data.frame(annotation),
                        params = params),
                   class = "GenomicRegion"))
}



# HELPERS

filter_by_gtf <- function(dt, list_of_args) {
  stopifnot(is.data.table(dt))


  dt <- as.data.table(dt)

  for (key in names(list_of_args)) {
    val <- list_of_args[[key]]

    # Ignore 1: If value not provided
    if (is.null(val) || all(is.na(val))) {
      message("Argument not given for ", key, " — skipping filter.")
      next
    }

    # Ignore 2: If col does not exist
    if (!(key %in% colnames(dt))) {
      message("Column ", key, " not found in data — skipping filter.")
      next
    }

    # Success: Perform filter
    dt <- dt[get(key) %in% val]
  }

  return(dt)
}




