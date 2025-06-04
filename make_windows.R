#!/usr/bin/Rscript
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(readr))


# Parse Arguments -----
upstream <- downstream <- out <- gtf <- bed <- biotype <- feature <- NULL
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    # Create the variable name (key) and get the value (val)
    key <- sub("^--", "", args[i])
    val <- if (i < length(args) && !startsWith(args[i + 1], "--")) args[i + 1] else TRUE
    
    # Bind the key to the value
    if (exists(key, inherits = FALSE)) {
      assign(key, val, envir = .GlobalEnv)
    } else { # if the argument name is not allowed
      stop(sprintf("Unknown argument: --%s", key))
    }
    
    i <- i + if (identical(val, TRUE)) 1 else 2
  } else {
    i <- i + 1
  }
}
if (any(sapply(c("upstream", "downstream", "out"), function(x) is.null(get(x))))) {
  stop("Missing required arguments: --upstream, --downstream, --out")
}
if (is.null(gtf) && is.null(bed)) {
  stop("Please provide either --gtf or --bed.")
}
if (!is.null(gtf) && !is.null(bed)) {
  stop("Mutually exclusive: --gtf and --bed.")
} else if (!is.null(bed)) {
  input <<- import(bed) %>%
    as.data.frame() %>%
    rename(gene_id = name) %>% as.data.frame()
} else {
  input <<- import(gtf) %>%
    as.data.frame()
}
# Cast numerical console arguments:
upstream <- as.numeric(upstream)
downstream <- as.numeric(downstream)

# Create the directory path in out if they don't exist
dir_path <- dirname(out)
if (dir_path != "." && !dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
# |__ args parsed ----- -----

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

# Create the .bed file 
input %>%
  
  # Filter feature (gene, variaton, ...) and biotype (protein_coding, lncRNA, ...)
                          # columns are quoted
  filter(using_if_given(., "gene_biotype", biotype)) %>%
  filter(using_if_given(., "type", feature)) %>%
  
  # Calculating upstream-downstream windows
  mutate(start = ifelse(strand == "+", start-1-upstream, start-1-downstream),
         end = ifelse(strand == "+", end+downstream, end+upstream),
         score = ifelse(is.na(score), ".", score)) %>%
  
  # remove regions where coordinates go over the edge
      # Alternative could be to mutate to zero if negative...
  filter(start >= 0 || end >= 0) %>%
  select(seqnames, start, end, gene_id, score, strand) %>%
  write_tsv(file.path(out), col_names = F)

