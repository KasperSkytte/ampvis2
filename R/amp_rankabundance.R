#' Rank abundance plot
#'
#' Generates a rank abundance curve (rank abundance vs cumulative read abundance).
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param order_group A vector to order the groups by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Genus"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_log (\emph{logical}) Log10-scale the x-axis. (\emph{default:} \code{FALSE})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @import ggplot2
#' @importFrom dplyr group_by mutate summarise
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table setkey
#'
#' @export
#'
#' @section Preserving relative abundances in a subset of larger data:
#' See \code{?\link{amp_subset_samples}} or the \href{https://madsalbertsen.github.io/ampvis2/articles/faq.html#preserving-relative-abundances-in-a-subset-of-larger-data}{ampvis2 FAQ}.
#'
#' @details A rank abundance curve is used to assess the biodiversity by plotting the ranked abundances of the OTUs (rank 1 is the most abundant, rank 2 the second and so on) versus the cumulative read abundance of the particular OTU. The rank abundances of the OTUs can be grouped by any taxonomic level by the \code{tax_aggregate} argument, the default is per OTU. When the samples are grouped by the \code{group_by} argument, the average of the group is used.
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Rank abundance plot
#' amp_rankabundance(AalborgWWTPs, group_by = "Plant")
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_rankabundance <- function(data,
                              group_by = "Sample",
                              order_group = NULL,
                              plot_log = FALSE,
                              tax_add = NULL,
                              tax_aggregate = "Genus",
                              tax_empty = "best",
                              tax_class = NULL,
                              normalise = TRUE,
                              detailed_output = FALSE) {

  ### Data must be in ampvis2 format
  is_ampvis2(data)

  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)

  # tax_add and tax_aggregate can't be the same
  if (!is.null(tax_aggregate) & !is.null(tax_add)) {
    if (tax_aggregate == tax_add) {
      stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
    }
  }

  # normalise counts
  if (isTRUE(normalise)) {
    data <- normaliseTo100(data)
  }

  ## Make a name variable that can be used instead of tax_aggregate to display multiple levels
  suppressWarnings(
    if (!is.null(tax_add)) {
      if (tax_add != tax_aggregate) {
        data$tax <- data.frame(data$tax, Display = apply(data$tax[, c(tax_add, tax_aggregate)], 1, paste, collapse = "; "))
      }
    } else {
      data$tax <- data.frame(data$tax, Display = data$tax[, tax_aggregate])
    }
  )

  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = data$tax[, "Display"], data$abund) %>%
    tidyr::gather(key = Sample, value = Abundance, -Display) %>%
    as.data.table()

  abund3 <- abund3[, "Abundance" := sum(Abundance), by = list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    as.data.frame()

  ## Add group information
  suppressWarnings(
    if (group_by != "Sample") {
      if (length(group_by) > 1) {
        grp <- data.frame(Sample = rownames(data$metadata), Group = apply(data$metadata[, group_by], 1, paste, collapse = " "))
      } else {
        grp <- data.frame(Sample = rownames(data$metadata), Group = data$metadata[, group_by])
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else {
      abund5 <- data.frame(abund3, Group = abund3$Sample)
    }
  )

  temp3 <- group_by(abund5, Display, Group) %>%
    summarise(Mean = mean(Abundance))

  TotalCounts <- temp3[with(temp3, order(-Mean)), ] %>%
    group_by(Group) %>%
    mutate(dummy = 1) %>%
    mutate(Cumsum = cumsum(Mean), Rank = cumsum(dummy)) %>%
    as.data.frame()

  if (!is.null(order_group)) {
    TotalCounts$Group <- factor(TotalCounts$Group, levels = rev(order_group))
  }

  if (group_by != "Sample") {
    p <- ggplot(TotalCounts, aes_string("Rank", "Cumsum", color = "Group"))
  }
  if (group_by == "Sample") {
    p <- ggplot(TotalCounts, aes_string("Rank", "Cumsum"))
  }
  p <- p +
    geom_line(size = 1) +
    ylim(0, 100) +
    xlab("Rank abundance") +
    theme_classic()

  if ((is.null(attributes(data)$normalised) | isFALSE(attributes(data)$normalised)) &
    isFALSE(normalise)) {
    p <- p + ylab("Cumulative read counts")
  } else {
    p <- p + ylab("Cumulative read abundance(%)")
  }

  if (plot_log == TRUE) {
    p <- p + scale_x_log10()
  }

  ## Define the output
  if (detailed_output) {
    return(list(plot = p, data = TotalCounts))
  } else if (!detailed_output) {
    return(p)
  }
}
