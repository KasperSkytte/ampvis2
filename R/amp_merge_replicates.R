#' @title Merge replicate samples
#'
#' @description Aggregates read counts in replicate samples by calculating the mean abundances of OTU's.
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param merge_var Variable in the sample metadata that defines the sample groups, see examples.
#' @param round If the read counts in \code{data$abund} are integers, any decimals present after merging will be rounded either \code{"up"} or \code{"down"}. Make sure this makes sense if the read counts have been normalised, as it may result in 0's, 1's, and 2's everywhere. (\emph{default:} \code{NULL})
#'
#' @return An object of class \code{ampvis2}
#'
#' @importFrom magrittr %>%
#' @importFrom data.table data.table dcast melt %chin%
#' @importFrom stats as.formula na.omit
#' @export
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#'
#' @examples
#' # load data
#' d <- amp_load(ampvis2::example_otutable, ampvis2::example_metadata)
#' d
#'
#' # add a grouping variable to the sample metadata defining the groups of sample replicates
#' d$metadata$group <- c("group1", "group1", "group2", "group2", "group3", "group4", NA, NA)
#' d$metadata
#'
#' # merge by "group" rounding up the resulting values
#' dmerged <- amp_mergereplicates(d,
#'   merge_var = "group",
#'   round = "up"
#' )
#' dmerged$metadata
#' dmerged
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
amp_merge_replicates <- function(data,
                                 merge_var,
                                 round = NULL) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  ### merge_var must be a length 1 string and present in data
  if (!is.character(merge_var) || length(merge_var) != 1) {
    stop("\"merge_var\" must be a character vector of length 1", call. = FALSE)
  }

  if (!any(colnames(data$metadata) == merge_var)) {
    stop(paste0("\"", merge_var, "\" was not found among metadata variable names"), call. = FALSE)
  }

  # warn when user rounds up or down if data is normalised
  if (isTRUE(attributes(data)$normalised) & any(tolower(round) %in% c("down", "up"))) {
    warning("The data has been normalised, rounding up or down likely does not make sense",
      call. = FALSE
    )
  }

  # find unique group names excluding NA and empty strings ("")
  data$metadata[[merge_var]] <- as.character(data$metadata[[merge_var]])
  groups <- data$metadata[which(data$metadata[, merge_var] != ""), merge_var]
  groups <- unique(na.omit(groups))

  # add a group column to the abundance table to define the groups
  tempabund <- data$metadata[, c(1, which(colnames(data$metadata) == merge_var))] %>%
    merge(as.data.frame(t(data$abund)), by.x = 1, by.y = 0, sort = FALSE) %>%
    {
      .[, 1] <- ifelse(.[, 2] %chin% groups, .[, 2], .[, 1])
      colnames(.)[1] <- merge_var
      .[, -2]
    } %>%
    data.table::data.table()

  # melt to long format, calculate mean per group and
  # taxon, and cast back to wide format
  newabund <- data.table::melt(tempabund, id.vars = merge_var)[
    ,
    .(value = mean(value)),
    by = c(merge_var, "variable")
  ] %>%
    data.table::dcast(as.formula(paste(merge_var, "variable", sep = "~")))

  # transpose the new abundance table, remove first column and round up or down
  out <- data
  out$abund <- newabund %>%
    as.data.frame() %>%
    {
      # first column to rownames
      rownames(.) <- .[[1]]
      .[[1]] <- NULL
      .
    } %>%
    t() %>%
    {
      if (is.null(round)) {
        .
      } else if (tolower(round) == "up") {
        ceiling(.)
      } else if (tolower(round) == "down") {
        floor(.)
      } else {
        .
      }
    } %>%
    as.data.frame()

  # create new metadata by removing duplicates and only keep
  # the first row in each new sample group. Move merge variable to first column
  # and use as new sample ID's
  out$metadata[[1]] <- tempabund[[1]]
  out$metadata <- out$metadata[
    ,
    colnames(out$metadata) != merge_var,
    drop = FALSE
  ]
  colnames(out$metadata)[1] <- merge_var
  out$metadata <- out$metadata[!duplicated(out$metadata[[1]]), ]
  rownames(out$metadata) <- out$metadata[[1]]

  # make sure order of samples are the same in abund and metadata
  out$abund <- out$abund[, rownames(out$metadata), drop = FALSE]
  return(out)
}

#' @rdname amp_merge_replicates
#' @export
amp_mergereplicates <- amp_merge_replicates
