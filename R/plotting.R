.makeRegions_plotting_env <- new.env(parent = emptyenv())
.makeRegions_plotting_env$plotting_lib <- NULL

source("../metagene/binGenome.lib.R")

#' Wrapper for sourcing plotting_lib.R
#'
#' @rdname plotting
#' @name plotting
#'
#' @param path Path to the plotting_lib.R
#'
#' @return TRUE on successful execution, FALSE otherwise.
#' @export
source_plotting_lib <- function(path) {
  if (!file.exists(path)) stop("Invalid binGenome.sh path: ", path)
  if (!endsWith(path, ".R")) stop("Invalid file. Should be an R script.")
  .makeRegions_plotting_env$plotting_lib <- normalizePath(path)
  source(path)
  invisible(TRUE)
}


#' @rdname plotting
#' @description
#' Basic plot drawing function wrapped from binGenome plotting library
#'
#' @param plot_dir Where to save PDF of the plot
#' @param file_name File identifier. "metagene_" gets added as prefix.
#' @param binMatrixList List of binMatrix instances.
#' @param pvThresholdsVec A vector of numeric between 0 and 1 to use
#' as p-value thresholds in defaultPvalueColorTransformer for Wilcoxon test.
#'
#' @export
draw_metagene <- function(plot_dir, file_name, binMatrixList, pvThresholdsVec = c(0.05, 10^-3, 10^-5), ...) {

  # Default arguments
  defaults <- list(
    palette = grDevices::rainbow(length(binMatrixList)),
    legendPos = "topleft",
    legendSpacingBetween = 20,
    title = "",
    normBylibSize = FALSE,
    normByShapeSum = TRUE,
    normByBinlength = FALSE,
    normByShapeMax = FALSE,
    labels = NULL,
    cexLegend = 2,
    cexLegendSymbol = 3,
    cexLab = 1.4,
    fixedLabelsEndTotalBins = 0,
    fixedLabelsEnd = NULL,
    fixedLabelsStartTotalBins = 0,
    fixedLabelsStart = NULL,
    aggregateFun = "mean",
    legendNCol = 1,
    legendXIntersp = 0.8,
    legendYIntersp = 1,
    minYA = 0
  )

  # Unite formal and dot arguments into list
  dots <- list(...)

  # Keep the relevant arguments for wrapped plot function and warn for unused
  all_arguments <- names(formals(plotGroup))
  extra_args <- setdiff(names(dots), all_arguments)
  if (length(extra_args) > 0) {
    warning("Ignoring unused arguments passed to draw_metagene: ", paste(extra_args, collapse = ", "))
  }
  dots <- dots[names(dots) %in% all_arguments]

  # Update the defaults by the user's valid arguments
  defaults[names(dots)] <- dots

  # Final arguments the call wrapped function
  final_args <- c(mget("binMatrixList", envir = environment()), defaults)

  pdf(paste0(plot_dir, "/metagene_", file_name, ".pdf"), width = 14, height = 7)

  if (length(binMatrixList) == 2) { # Wilcox test if there are two conditions
    final_args$performWilcoxTest <- TRUE
    final_args$wilcoxTestPVTransformFUN <- function(x) defaultPvalueColorTransformer(x, pvThresholds = pvThresholdsVec)
  } else {
    final_args$performWilcoxTest <- FALSE
  }

  message("Calling plotGroup with ", paste(str(final_args), collapse = ", "))
  do.call(plotGroup, final_args)

  dev.off()
}



#' Pairwise and all
#'
#' Plots all pairwise comparisons and one with all conditions.
#' Exploits multithreading from package parallel using mclapply on pairs generated.
#' binMatrix objects in the coverage folder in run_dir are aggregated
#' using aggregate_FUN (default: mean) and saved as rds objects in a directory
#' called 'aggregated_binmatrices' as a cache for further calls of the same function
#' with same run_dir, annotation and aggregate function.
#'
#' @param plot_dir Where to save PDF outputs of plotting
#' @param experiment_dir The run folder where windows and binning took place.
#' The folder should have a 'config.json' file where the parameters for make_windows (optional)
#' and bin_genome (required) are given.
#' @param annotation_name Indicates the annotation used for binning. Appends to the file name.
#' @param condition_mapper A data frame with columns "grep_name", "color", "display_name"
#' to match conditions with line colors and names to display in the legend.
#' @param threads Number of threads to use in pairwise plotting.
#' Usage makes sense if multiple pairs are plotted, e.g. nrow(condition_mapper) > 1,
#' and cores are available, e.g. if detectCores() >= threads or
#' upstream threading has not used them up. Since this utilizes mclapply from
#' package parallel, this option should be used careful with nested parallelization.
#'
#' @returns
#' @export
#'
#' @examples
plot_metagene_experiment <- function(plot_dir, run_dir, annotation_name, condition_mapper,
                            aggregate_FUN = mean, threads = 1) {
  # Check arguments
  if(!dir.exists(run_dir)) stop("Run dir does not exist")
  if(!file.exists(file.path(run_dir, 'config.json'))) stop("Run dir does not have a config file.")
  absent_columns_in_mapper <- setdiff(c("grep_name", "color", "display_name"), names(condition_mapper))
  if(length(absent_columns_in_mapper) > 0)
    stop("condition_mapper data frame does not have the columns: ", paste(absent_columns_in_mapper, collapse = ", "))
  if(!dir.exists(plot_dir)) dir.create(plot_dir)

  # Rely on config to understand run structure
  config <- jsonlite::read_json(file.path(run_dir, "config.json"), simplifyVector = TRUE)
  cov_folder <- file.path(config$bin_genome$output_dir)
  condition_greps <- condition_mapper[["grep_name"]]

  message("Started processing ", run_dir, " with annotation:", annotation_name)

  # === label settings ===
  firstLabel <- "0%"
  endLabel   <- "100%"

  # include 0 and 100 if plotting only gene body
  # upstream labels
  if ("upstream" %in% config$make_windows) {
    firstLabel   <- ""
    upstream_kb  <- as.numeric(config$make_windows$upstream)
    up_by        <- upstream_kb / 5
    fixedLabelsStart <- c(
      paste0(seq(-upstream_kb, -1, by = up_by), "kb"),
      "start"
    )
  }
  # downstream labels
  if ("downstream" %in% config$make_windows) {
    endLabel     <- ""
    downstream_kb <- as.numeric(config$make_windows$downstream)
    down_by      <- downstream_kb / 5
    fixedLabelsEnd <- c(
      "end",
      paste0(seq(0, downstream_kb, by = down_by), "kb")
    )
  }

  # body labels
  bodyLabels <- c(firstLabel, paste0(seq(20, 80, by = 20), "%"), endLabel)

  # fixed bins
  fixedLabelsEndTotalBins <- ifelse(
    length(config$bin_genome$fixedBinSizeDownstream) == 0,
    0L,
    as.integer(strsplit(as.character(config$bin_genome$fixedBinSizeDownstream), ":")[[1]][1])
  )
  fixedLabelsStartTotalBins <- ifelse(
    length(config$bin_genome$fixedBinSizeDownstream) == 0,
    0L,
    as.integer(strsplit(as.character(config$bin_genome$fixedBinSizeUpstream), ":")[[1]][1])
  )


  # Create pairs
  pairs <- combn(condition_greps, 2) # iterate through cols for pairwise comparisons

  # === binMatrix aggregation ===
  # Check if grouped binMatrices exist as Rds objects
  aggregated_matrices <- list()
  rds_dir <- "aggregated_binmatrices"
  aggr_func <- deparse(substitute(aggregate_FUN))
  rds_object <- file.path(rds_dir,
                          paste0(paste0(strsplit(run_dir, "/")[[1]], collapse = "_"),
                                 ".Rds") #"_", _aggr_func,
  )
  if (!dir.exists(rds_dir) & !file.exists(rds_object)) {
    message("No RDS object found, reading and grouping and saving as RDS for next time!..")
    # read bin tables from the coverage folder
    bin_tables <- list()
    for (cond_index in seq_along(condition_greps)) {
      cond_name <- condition_greps[[cond_index]]

      bin_tables[[cond_name]] <- getBinTable(folder = cov_folder,
                                             condition = cond_name,
                                             ignore = ".norm.coverage.csv",
                                             gsub_name = "aggregate_gene")
    }

    # Assign names for named access: bin_tables[[display_name]]
    names(bin_tables) <- condition_mapper$display_name
    # summarize each condition using median of reps
    grouped <- lapply(bin_tables, function(condition) {
      groupMatrices(condition, median)
      }
    )
    dir.create(rds_dir)
    saveRDS(grouped, file = rds_object)
  } else {
    message("Reading from the found RDS object: ", rds_object)
    grouped <- readRDS(rds_object)
  }

  # plot pairs
  message("Pairwise plotting in threads: ", threads)
  mclapply(1:ncol(pairs), function(comparison_index) {
    this <- pairs[1, comparison_index]
    that <- pairs[2, comparison_index]

    # match preserves order of this - that
    color_idx <- match(c(this, that), condition_mapper$grep_name)
    colors <- condition_mapper$color[color_idx]

    name_idx <- match(c(this, that), condition_mapper$grep_name)
    display_names <- condition_mapper$display_name[name_idx]

    message("Drawing plot for ", paste(display_names, collapse = ", "),
                  " with color(s) ", paste(colors, collapse = ", "), ".")

    pairwise_comparison <- grouped[display_names]

    draw_metagene(
      plot_dir,
      paste0(annotation_name, "_", paste0(gsub(" ", "_", display_names), collapse = "_vs_")),
      pairwise_comparison,
      labels = bodyLabels,
      ylab = "Coverage",
      legendPos = "top",
      title = paste(display_names, collapse = " vs "),
      palette = colors,
      fixedLabelsStartTotalBins = fixedLabelsStartTotalBins,
      fixedLabelsEndTotalBins = fixedLabelsEndTotalBins
    )

    return(paste(display_names, collapse = "_vs_"))
  }, mc.cores = threads, mc.preschedule = TRUE, mc.cleanup = TRUE)

  # Plot all genes:
  draw_metagene(plot_dir, paste0("all_", annotation_name), grouped,
            labels = bodyLabels,
            ylab = "Coverage",
            legendPos = "top",
            palette = condition_mapper$color,
            fixedLabelsStartTotalBins = fixedLabelsStartTotalBins,
            fixedLabelsEndTotalBins = fixedLabelsEndTotalBins)
}
