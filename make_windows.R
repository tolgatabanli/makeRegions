#!/usr/bin/Rscript
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(readr))


# Parse Arguments -----
biotype <- NA
upstream <- downstream <- out <- gtf <- bed <- NULL
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    key <- sub("^--", "", args[i])
    val <- if (i < length(args) && !startsWith(args[i + 1], "--")) args[i + 1] else TRUE
    
    if (exists(key, inherits = FALSE)) {
      assign(key, val, envir = .GlobalEnv)
    } else {
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
upstream <- as.numeric(upstream)
downstream <- as.numeric(downstream)

dir_path <- dirname(out)
if (dir_path != "." && !dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
# ----- -----

# Create the .bed file 
input %>%
  filter(ifelse(is.na(biotype) | !("gene_biotype" %in% colnames(.)),
                TRUE,
                gene_biotype == biotype)) %>%
  mutate(start = ifelse(strand == "+", start-1-upstream, start-1-downstream),
         end = ifelse(strand == "+", end+downstream, end+upstream),
         score = ifelse(is.na(score), ".", score)) %>%
  select(seqnames, start, end, gene_id, score, strand) %>%
  write_tsv(file.path(out), col_names = F)

