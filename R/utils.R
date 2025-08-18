# Manage namespaces and utility functions

#' @importFrom magrittr %>%
#' @importFrom data.table := fifelse
NULL

.makeRegions_env <- new.env(parent = emptyenv())

.makeRegions_env$configFile <- NULL
.makeRegions_env$logFile <- NULL

#' Set config and log file paths for makeRegions
#'
#' These functions set and create the files for config (.json) and log (.log)
#' files used by the \code{makeRegions} process.
#'
#' \itemize{
#'   \item \code{setConfig()} sets the file path for the configuration file.
#'   \item \code{setLog()} sets the file path for the log file.
#' }
#'
#' @param file Character string specifying the file path. The containing
#'   directory will be created if it does not exist. Use \code{NULL} to unset.
#'
#' @returns No return value, called for side effects.
#' @export
#' @rdname makeRegionsFiles
#'
#' @examples
#' \dontrun{
#' setConfig("~/project/config.json")
#' setLog("~/project/console.log")
#' }
setConfig <- function(file) {
  if (!is.null(file)) {
    stopifnot(is.character(file), length(file) == 1)
    dir_path <- dirname(file)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    if(!file.exists(file)) {
      file.create(file)
    }
    writeLines("{}", file) # initiate json
  }
  .makeRegions_env$configFile <- file
}

#' @export
#' @rdname makeRegionsFiles
setLog <- function(file) {
  if (!is.null(file)) {
    stopifnot(is.character(file), length(file) == 1)
    dir_path <- dirname(file)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    if(!file.exists(file)) {
      file.create(file)
    }
  }
  .makeRegions_env$logFile <- file
}

write_config <- function(caller_name, params) {
  # write the parameters of the function into the config file
  if (!is.null(.makeRegions_env$configFile)) {
    config <- jsonlite::read_json(.makeRegions_env$configFile)
    config[[caller_name]] <- params
    jsonlite::write_json(config, path = .makeRegions_env$configFile,
                         pretty = TRUE, auto_unbox = TRUE)
  }
}

# caller_name: function name
# args: named list of parameter and values
start_log <- function(caller_name, args) {
  if (!is.null(.makeRegions_env$logFile)) {
    log_con <- file(.makeRegions_env$logFile, open = "at")

    sink(log_con, type = "output", append = TRUE)
    sink(log_con, type = "message", append = TRUE)

    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(paste0("\n[", timestamp, "]"))

    args <- args[!vapply(args, is.null, FUN.VALUE = logical(1))]
    args <- paste(
      vapply(names(args), function(n) {
        val <- args[[n]]
        # Keep it short for large objects
        if (is.character(val) || is.numeric(val) || is.logical(val)) {
          paste0(n, " = ", toString(val))
          } else {
            paste0(n, " = <", class(val)[1], ">")
            }
        }, FUN.VALUE = character(1)),
      collapse = ", ")
    message(paste0(caller_name, "(", args, ")"))

    return(function() {
      sink(type = "output")
      sink(type = "message")
      close(log_con)
    })
  }
}


