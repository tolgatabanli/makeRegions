utils::globalVariables(c("gtf", "bed", "name", "path_to_out", "start", "end", "strand", "seqnames", "score", "input", ":="))

#' Create genomic regions with upstream-downstream sizes and optional filters.
#'
#' @param input_file GTF/GFF or BED file
#' @param upstream Upstream size
#' @param downstream Downstream size
#' @param path_to_output (Optional) File path to write output. Creates dirs if they do not exist
#' @param position (Optional) Specify 'start' or 'end' to build the window around.
#' If not specified, the start and end coordinates are extended with upstream and downstream sizes, correspondingly
#' @param ... (Optional) GTF arguments to filter. Only if GTF file is given.
#' For example: type = c("gene", "transcript", "exon"...),
#' gene_biotype = c("lncRNA", "pseudogene"...) etc.
#' A filter argument with a vector of length 2 or more chooses rows that satisfy ANY one of the elements.
#'
#' @returns A data.frame in .bed format with specified window and filter criteria.
#' Columns are "seqnames", "start", "end", "name", "score", "strand".
#' 'seqnames' and 'strand' are factors.
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
  params <- as.list(environment())
  score <- NA
  gtf_filters <- list(...)

  # Argument checks
  if (missing(input_file) || missing(upstream) || missing(downstream)) {
    stop("Missing required arguments: upstream, downstream")
  }

  is_gtf <- FALSE # use for ID column creation
  if(endsWith(input_file, ".gtf") | endsWith(input_file, ".gff")) {
    is_gtf <- TRUE
    imported_file <- rtracklayer::import(input_file) %>% # seqnames, start, end, width, strand, source, type, score, phase, gene_id ...
      as.data.frame()
  } else if (endsWith(input_file, ".bed")) {
    imported_file <- rtracklayer::import(input_file) %>% # seqnames, start, end, width, strand, name, score
      as.data.frame()
    if(length(gtf_filters) > 0) {
      stop("GTF filters cannot be used with input file of format BED.")
    }
  } else {
    stop(paste0("Provided input file (", input_file, ") is of wrong format. Should: .gtf, .gff or .bed"))
  }

  # optional argument check
  if(!is.null(position) && !(position %in% c("start", "end"))) {
    stop("Position argument should be either 'start' or 'end' !")
  }

  # Create the directory path in out if they don't exist
  if (!is.null(path_to_output)) {
    dir_path <- dirname(path_to_output)
    if (dir_path != "." && !dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  # ================ Config and Logs  ================
  write_config("make_windows", params)
  revert_sink <- start_log("make_windows", params) # starts the log and returns a function to revert the sink on exit
  if (!is.null(revert_sink)) on.exit(revert_sink(), add = T)


  #  ========================== Creating BED ==========================

  annotation <- data.table::as.data.table(imported_file)

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
    start = fifelse(strand == "+", start - 1 - upstream, start - 1 - downstream),
    end   = fifelse(strand == "+", end + downstream, end + upstream)
  )]

  # if score not found, fill with dot in .bed annotation
  annotation[, score := as.character(score)]
  annotation[is.na(score), score := "."]

  # Remove regions where coordinates go over the edge with a message
  negative_detected <- annotation[, start < 0 | end < 0]
  if (any(negative_detected, na.rm = TRUE)) {
    message(paste("Negative coordinates filtered out!", sum(negative_detected, na.rm = TRUE), "rows eliminated"))
  }

  annotation <- annotation[!negative_detected]

  # Create name column based on identifiers if GTF is given
  if (is_gtf) {
    annotation[, name := NA_character_]
    annotation[type == "gene" | type == "aggregate_gene", name := gene_id]
    annotation[type == "transcript", name := paste(gene_id, transcript_id, sep = "-")]
    annotation[type == "exonic_part", name := paste(gene_id, exonic_part_number, sep = "-")]
    annotation[type == "exon",       name := paste(gene_id, transcript_id, exon_number, sep = "-")]

    annotation <- annotation[, .(seqnames, start, end, name, score, strand)]
  }
  annotation <- annotation[, .(seqnames, start, end, name, score, strand)]



  # row number for tracking changes and logging
  message(paste0("Input had ", nrow(imported_file),
                " rows, output has ", nrow(annotation), " rows (",
                data.table::uniqueN(annotation), " unique)."))

  if(!is.null(path_to_output)) {
    annotation %>%
      readr::write_tsv(file.path(path_to_output), col_names = F)
    message(paste("File has been written:", path_to_output))
  }

  invisible(annotation)
}



# HELPERS

filter_by_gtf <- function(dt, list_of_args) {
  stopifnot(is.data.table(dt))


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




