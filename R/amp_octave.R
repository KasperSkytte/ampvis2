#' @title Octave plot
#'
#' @description Generates an octave plot to assess diversity. An octave plot is a histogram of the number of OTU's observed by bins of read counts, where the bin ranges increase exponentially, see details.
#'
#' @usage amp_octave(data, facet_by = "")
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param facet_by Split the plot into subplots based on a categorical/discrete variable in the samples metadata. (\emph{default:} \code{NULL})
#' @param scales If \code{facet_by} is set, should the axis scales of each subplot be fixed (\code{fixed}), free (\code{"free"}), or free in one dimension (\code{"free_x"} or \code{"free_y"})? (\emph{default:} \code{"fixed"})
#' @param num_threads Maximum number of distinct groups as defined by \code{facet_by} to process simultaneously using multicore processors. Only used if \code{facet_by} is set and there are more than one distinct group in the given variable(s). Defaults is the number of available cores minus 1.
#'
#' @details The \eqn{n}th bin in the histogram has the range \eqn{r(n)=2^n...2^{n+1}-1}.
#' The height of the bars then reflect the number of unique OTU's with read counts in each bin.
#' By judging the distribution one can assess whether the samples have been sequenced deeply
#' enough. A full symmetrical bell-shaped distribution with the left part far from the y-axis is the ideal.
#' A high amount of OTU's with a low amount of reads indicates noise, chimeras, and even cross talk.
#'
#' @return A ggplot2 object
#' @importFrom dplyr as_tibble tibble
#' @importFrom foreach foreach %do%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom plyr ddply
#' @import ggplot2
#'
#' @export
#'
#' @seealso
#' \code{\link{amp_load}}, \code{\link{amp_rarecurve}}
#'
#' @examples
#' # Load example data
#' data("MiDAS")
#'
#' # Generate an octave plot of all samples
#' amp_octave(MiDAS)
#'
#' # Generate an octave plot for each group as
#' # defined by the metadata variable "Year", set the
#' # y-axes free and process 4 groups simultaneously using
#' # multicore processing
#' amp_octave(MiDAS,
#'   facet_by = "Year",
#'   scales = "free_y",
#'   num_threads = 4
#' )
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
amp_octave <- function(data,
                       facet_by = NULL,
                       scales = "fixed",
                       num_threads = parallel::detectCores() - 1L) {
  abund <- data$abund
  # check if samples in metadata and abund match, and that their order is the same
  if (!identical(colnames(abund), data$metadata[[1]])) {
    if (!all(colnames(abund) %in% data$metadata[[1]]) | !all(data$metadata[[1]] %in% colnames(abund))) {
      stop("The samples in metadata do not match those in the read counts table!", call. = FALSE)
    }
    warning("The order of samples in the sample metadata ($metadata) and in the read
            counts table ($abund) are not the same, reordering read counts table...", call. = FALSE)
    # reorder abund based on metadata
    abund <- abund[, data$metadata[[1]], drop = FALSE]
  }

  # transpose abund and merge with one or more metadata variable(s) as defined by facet_by
  OTUsums <- dplyr::as_tibble(t(abund))
  if (!is.null(facet_by)) {
    # check if metadata variable exists (both by name or position)
    if ((is.character(facet_by) & !all(facet_by %in% colnames(data$metadata))) |
      (is.numeric(facet_by) & !all(facet_by >= 1 & facet_by <= ncol(data$metadata)))) {
      stop("One or more of the metadata variables(s) provided 
           with facet_by was not found in the sample metadata!", call. = FALSE)
    }

    OTUsums$group <- apply(data$metadata[, facet_by, drop = FALSE], 1, paste, collapse = ", ")
    # multiprocessing is done by group, so only relevant if there is more than one distinct
    # group in the variable(s) set with facet_by
    if (length(unique(OTUsums$group)) > 1L & is.numeric(num_threads) & num_threads > 1L) {
      doParallel::registerDoParallel(num_threads, num_threads)
      .parallel <- TRUE
    }
  } else {
    # all samples at once if facet_by is not set
    OTUsums$group <- "All samples"
    .parallel <- FALSE
  }

  # calculate OTU sums for each group and bin them
  # Several different implementations to calculate OTU sums per group was tested
  # dplyr with group_by(group) and summarise_all(sum) was extremely slow,
  # and data.table was even slower, with fx OTUsums[,.(nOTUs = colSums(.SD)), by = group]
  # or OTUsums[,.(nOTUs = lapply(.SD, sum)), by = group] and variations thereof.
  # Using base::colSums and a sequential for loop to generate the bin sums was the fastest
  # ddply calculates group in parallel, dont use %dopar% in the foreach loop, its slower,
  # its only 15-20 bins anyways. Cant use apply family of functions as the bin
  # ranges have to be calculated based on the next bin. Generates a long format table to be
  # able to facet with ggplot2
  binSums <- plyr::ddply(
    .data = OTUsums,
    .variables = "group",
    .parallel = .parallel,
    .fun = function(x) {
      OTUsums <- colSums(x[, -which(colnames(x) == "group")])
      bins <- 2^(0:ceiling(log2(max(OTUsums)))) # generate bin sizes based on data
      binSums <- foreach::foreach(
        i = seq_along(bins),
        .combine = "rbind",
        .inorder = TRUE
      ) %do% {
        # sum the number of OTUs for each bin with range 2^n...2^(n+1)-1,
        # meaning including lower bound, but not upper bound as it is
        # the start of the next bin
        dplyr::tibble(
          bin = bins[i],
          nOTUs = sum(OTUsums >= bins[i] &
            OTUsums < bins[i + 1])
        )
      }
      binSums$bin <- as.factor(binSums$bin)
      binSums
    }
  )
  if (isTRUE(.parallel)) {
    doParallel::stopImplicitCluster()
  }

  # gogo plot, facet if facet_by is not NULL
  plot <- ggplot(binSums, aes(bin, nOTUs)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Minimum reads") +
    ylab("OTU's") + {
      if (!is.null(facet_by)) {
        facet_wrap(~group, scales = scales)
      }
    }
  return(plot)
}
