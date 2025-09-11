#!/home/proj/software/R/R-4.2.2.mkl/bin/Rscript

source("../metagene/plotting_lib.R")

#' Wrapper for basic plot drawing from binGenome plotting library
#'
#' @param plot_folder Where to save PDF of the plot
#' @param file_name File identifier. "metagene_" gets added as prefix.
#' @param binMatrixList List of binMatrix instances.
#'
#' @export
#'
#' @examples
draw_plot <- function(plot_folder, file_name, binMatrixList ...) {
  # Default arguments
  defaults <- list(
    palette = grDevices::rainbow(length(binMatrixList)),
    legendPos = "topleft",
    legendSpacingBetween = 20,
    title = "",
    pvThresholdsVec = c(0.05, 10^-3, 10^-5),
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

  # unite formal and dot arguments into list
  dots <- list(...)
  args <- c(mget(names(formals()), envir = environment()), dots)

  pdf(paste0(plot_folder, "/metagene_", file_name, ".pdf"), width = 14, height = 7)


  # Keep the relevant arguments for wrapped plot function and warn for unused
  required_args <- names(formals(plotGroup))
  extra_args <- setdiff(names(args), required_args)
  if (length(extra_args) > 0) {
    warning("Ignoring unused arguments passed to draw_plot: ", paste(extra_args, collapse = ", "))
  }
  args <- args[names(args) %in% required_args]


  if (length(binMatrixList) == 2) { # Wilcox test if there are two conditions
    args$performWilcoxTest <- TRUE
    args$wilcoxTestPVTransformFUN <- function(x) defaultPvalueColorTransformer(x, pvThresholds = pvThresholdsVec)
  } else {
    args$performWilcoxTest <- FALSE
  }

  do.call(plotGroup, args)

  dev.off()
}


# Plots pairwise and all comparisons
#
# @param plot_dir Folder to save the plots as .pdf
# @param experiment_folder The folder where windows and binning took place (~ "igv/whole_gene" etc.)
# @param condition_mapper A data frame with columns "grep_name", "color", "display_name"
plot_experiment <- function(plot_dir, experiment_folder, region, condition_mapper, remaining_cores) {
  # Rely on config to understand run structure
  config <- jsonlite::read_json(file.path(experiment_folder, "config.json"), simplifyVector = TRUE)
  cov_folder <- file.path(config$bin_genome$output_dir)
  condition_greps <- condition_mapper[["grep_name"]]

  message(paste("Started processing", experiment_folder, "with region:", region))

  # === label settings ===
  firstLabel <- "0%"
  endLabel   <- "100%"

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
    is.na(config$bin_genome$fixedBinSizeDownstream),
    0L,
    as.integer(strsplit(as.character(config$bin_genome$fixedBinSizeDownstream), ":")[[1]][1])
  )
  fixedLabelsStartTotalBins <- ifelse(
    is.na(config$bin_genome$fixedBinSizeUpstream),
    0L,
    as.integer(strsplit(as.character(config$bin_genome$fixedBinSizeUpstream), ":")[[1]][1])
  )

  # Create pairs of comparisons
  pairs <- combn(condition_greps, 2) # should iterate through cols for pairwise comparisons

  # === binMatrix aggregation ===
  # Check if grouped binMatrices exist as Rds objects
  aggregated_matrices <- list()
  rds_dir <- "aggregated_binmatrices"
  rds_object <- file.path(rds_dir,
                          paste0(paste0(strsplit(experiment_folder, "/")[[1]], collapse = "_"), ".Rds")
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
    message("Reading from the found RDS object!..")
    grouped <- readRDS(rds_object)
  }

  # plot pairs
  message(paste0("Future apply working with cores:", remaining_cores))
  mclapply(1:ncol(pairs), function(comparison_index) {
    this <- pairs[1, comparison_index]
    that <- pairs[2, comparison_index]

    # match preserves order of this - that
    color_idx <- match(c(this, that), condition_mapper$grep_name)
    colors <- condition_mapper$color[color_idx]

    name_idx <- match(c(this, that), condition_mapper$grep_name)
    display_names <- condition_mapper$display_name[name_idx]

    message(paste("Drawing plot for", paste(display_names, collapse = ", "),
                  "with colors", paste(colors, collapse = ", "), "respectively."))

    pairwise_comparison <- grouped[display_names]

    draw_plot(
      plot_dir,
      paste0(region, "_", paste0(gsub(" ", "_", display_names), collapse = "_vs_")),
      pairwise_comparison,
      labels = bodyLabels,
      ylab = "Coverage",
      legendPos = "top",
      title = paste(display_names, collapse = " vs "),
      palette = colors
    )

    return(paste(display_names, collapse = "_vs_"))
  }, mc.cores = remaining_cores, mc.preschedule = TRUE, mc.cleanup = TRUE)

  # Plot all genes:
  draw_plot(plot_dir, paste0("all_", region), grouped,
            labels = bodyLabels,
            ylab = "Coverage",
            legendPos = "top",
            palette = condition_mapper$color)
}
