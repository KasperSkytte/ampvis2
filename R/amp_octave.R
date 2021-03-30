#' @title Octave plot
#'
#' @description Generates an octave plot to assess alpha diversity. An octave plot is a histogram of the number of taxa observed by bins of read counts, where the bin ranges increase exponentially, see details.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_aggregate Aggregate (sum) OTU's to a specific taxonomic level initially. OTU's that have not been assigned at the chosen level will be removed with a message. (\emph{default:} \code{"OTU"})
#' @param group_by Group the samples based on a categorical/discrete variable in the metadata. It is recommended to look at samples individually. Can be a character vector with variable names as-is or a numerical vector with variable positions in the metadata of any length. Set to \code{NULL} for grouping all samples together (not recommended). (\emph{default:} \code{1})
#' @param scales If \code{group_by} is set, should the axis scales of each subplot be fixed (\code{fixed}), free (\code{"free"}), or free in one dimension (\code{"free_x"} or \code{"free_y"})? (\emph{default:} \code{"fixed"})
#' @param num_threads Maximum number of distinct groups as defined by \code{group_by} to process simultaneously using multicore processors. Only used if \code{group_by} is set and there are more than one distinct group in the given variable(s). Default is the number of available cores minus 1.
#'
#' @details The \eqn{n}th bin in the histogram has the range \eqn{r(n)=2^n...2^{n+1}-1}.
#' The height of the bars then reflect the number of unique taxa with read counts in each bin.
#' By judging the distribution one can assess whether the samples have been sequenced deeply
#' enough at the chosen taxonomic level. A full symmetrical bell-shaped distribution
#' with the left part far from the y-axis is the ideal.
#' A high amount of OTU's with a low amount of reads indicates noise, chimeras, and even cross talk.
#'
#' Aggregating OTU's using \code{tax_aggregate} is useful to assess whether the samples
#' have been sequenced deep enough to capture the full diversity at the given level, but
#' ONLY applies to OTU's that have assigned taxonomy at the given level.
#'
#' It is recommended to look at samples individually as grouping samples will almost always look ideal.
#' It is better to identify "bad" samples individually and remove them.
#'
#' @return A ggplot2 object
#' @importFrom data.table getDTthreads setDTthreads melt
#' @importFrom parallel detectCores
#' @import ggplot2
#'
#' @export
#'
#' @seealso
#' \code{\link{amp_load}}, \code{\link{amp_rarecurve}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Subset data
#' ds <- amp_subset_samples(AalborgWWTPs, Year %in% 2014)
#'
#' # Generate an octave plot of all samples at Genus level. Adjust num_threads to
#' # process multiple groups simultaneously using multicore processing
#' amp_octave(ds,
#'   group_by = "SampleID",
#'   tax_aggregate = "Genus",
#'   scales = "free_y",
#'   num_threads = 1
#' )
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
amp_octave <- function(data,
                       tax_aggregate = "OTU",
                       group_by = 1L,
                       scales = "fixed",
                       num_threads = parallel::detectCores() - 2L) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # function only works with counts, not relative abundances
  if (!abundAreCounts(data)) {
    stop("amp_octave can only be used with read counts, not relative abundances", call. = FALSE)
  }

  # check if samples in metadata and abund match, and that their order is the same
  if (!identical(colnames(data$abund), data$metadata[[1]])) {
    if (!all(colnames(data$abund) %in% data$metadata[[1]]) | !all(data$metadata[[1]] %in% colnames(data$abund))) {
      stop("The samples in metadata do not match those in the read counts table!", call. = FALSE)
    }
    warning("The order of samples in the sample metadata ($metadata) and in the read counts table ($abund) are not the same, reordering read counts table...", call. = FALSE)
    # reorder abund based on metadata
    data$abund <- data$abund[, data$metadata[[1]], drop = FALSE]
  }

  # check if metadata variable exists (both by name or position)
  if (!is.null(group_by)) {
    if ((is.character(group_by) & !all(group_by %in% colnames(data$metadata))) |
      (is.numeric(group_by) & !all(group_by >= 1 & group_by <= ncol(data$metadata)))) {
      stop("One or more of the metadata variables(s) provided with group_by was not found in the sample metadata!", call. = FALSE)
    }
  }

  # max number threads to use, data.table only
  DTthreads <- data.table::getDTthreads()
  if (is.numeric(num_threads)) {
    if (num_threads > 0L) {
      data.table::setDTthreads(threads = num_threads)
    }
  }

  # Aggregate OTUs to a specific taxonomic level
  lowestTaxLevel <- getLowestTaxLvl(data$tax, tax_aggregate)
  abundAggr <- aggregate_abund(
    abund = data$abund,
    tax = data$tax,
    tax_aggregate = lowestTaxLevel,
    tax_add = NULL,
    calcSums = FALSE,
    format = "long"
  )

  # Calculate the sum of each taxa by each group in group_by if set
  if (!is.null(group_by)) {
    # Merge with the variable(s) defined by group_by
    # (note: this is ~58 times faster than data.table::merge() by sample)
    abundAggr[,
      group := apply(
        data$metadata[which(data$metadata[[1]] == Sample), group_by, drop = FALSE],
        1,
        paste,
        collapse = ", "
      ),
      by = Sample
    ]
  } else {
    # all samples at once if group_by is not set
    abundAggr[, group := "All samples"]
  }

  # Count OTU's per group and tax level
  taxSums <- abundAggr[, .(taxSums = sum(Abundance)), by = .(Display, group)]

  # Count the number of distinct axa per bin
  # not at all concise, but it's by far the fastest
  binSums <- data.table::melt(
    id.vars = "group",
    variable.name = "bin",
    value.name = "nTaxa",
    variable.factor = TRUE,
    taxSums[, .(
      `1` = sum(taxSums == 1),
      `2` = sum(taxSums >= 2 & taxSums < 4),
      `4` = sum(taxSums >= 4 & taxSums < 8),
      `8` = sum(taxSums >= 8 & taxSums < 16),
      `16` = sum(taxSums >= 16 & taxSums < 32),
      `32` = sum(taxSums >= 32 & taxSums < 64),
      `64` = sum(taxSums >= 64 & taxSums < 128),
      `128` = sum(taxSums >= 128 & taxSums < 256),
      `256` = sum(taxSums >= 256 & taxSums < 512),
      `512` = sum(taxSums >= 512 & taxSums < 1024),
      `1024` = sum(taxSums >= 1024 & taxSums < 2048),
      `2048` = sum(taxSums >= 2048 & taxSums < 4096),
      `4096` = sum(taxSums >= 4096 & taxSums < 8192),
      `8192` = sum(taxSums >= 8192 & taxSums < 16384),
      `16384` = sum(taxSums >= 16384 & taxSums < 32768),
      `32768` = sum(taxSums >= 32768 & taxSums < 65536),
      `65536` = sum(taxSums >= 65536 & taxSums < 131072),
      `131072` = sum(taxSums >= 131072 & taxSums < 262144),
      `262144` = sum(taxSums >= 262144 & taxSums < 524288),
      `524288` = sum(taxSums >= 524288 & taxSums < 1048576)
    ),
    by = group
    ]
  )

  # gogo plot, facet if group_by is not NULL
  plot <- ggplot(binSums, aes(bin, nTaxa)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Minimum reads") +
    {
      if (tax_aggregate == "OTU") {
        ylab("Number of distinct OTU's")
      } else {
        ylab(paste0("Number of distinct taxa (", lowestTaxLevel, " level)"))
      }
    } +
    {
      if (!is.null(group_by)) {
        facet_wrap(~group, scales = scales)
      }
    }

  # restore previous DTthreads setting
  data.table::setDTthreads(threads = DTthreads)

  return(plot)
}
