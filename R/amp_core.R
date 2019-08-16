#' Core community analysis
#'
#' Generates a core community plot.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the samples by a variable in the metadata. (\emph{default:} \code{"Sample"})
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"OTU"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param abund_thrh Threshold in percent for defining "abundant"/"core" taxa. (\emph{default:} \code{0.1})
#' @param plotly (\emph{logical}) Returns an interactive plot instead. (\emph{default:} \code{FALSE})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @import ggplot2
#' @importFrom dplyr filter group_by mutate summarise
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table setkey data.table
#' @importFrom plotly ggplotly layout
#'
#' @export
#'
#' @details Plots the number of samples a given OTU is observed in (first axis) against the number of samples it is abundant in (second axis), as defined by the threshold. By setting the \code{group_by} argument to a variable in the metadata, then the axes will show their frequency in this group instead of per sample.
#'
#' The OTU points can be aggregated by the taxonomy by setting fx \code{tax_aggregate = "Phylum"}. To see the corresponding taxonomy of the OTUs use \code{plotly = TRUE} for an interactive plot.
#'
#' @section Preserving relative abundances in a subset of larger data:
#' See \code{?\link{amp_subset_samples}} or the \href{https://madsalbertsen.github.io/ampvis2/articles/faq.html#preserving-relative-abundances-in-a-subset-of-larger-data}{ampvis2 FAQ}.
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @references
#'   Saunders, Aaron M; Albertsen, Mads; Vollertsen, Jes; Nielsen, Per H (2016): The activated sludge ecosystem contains a core community of abundant organisms, ISME Journal 10, 11-20. \url{https://doi.org/10.1038/ismej.2015.117}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Core plot
#' amp_core(AalborgWWTPs)
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_core <- function(data,
                     group_by = "Sample",
                     abund_thrh = 0.1,
                     tax_aggregate = "OTU",
                     tax_class = NULL,
                     tax_empty = "best",
                     plotly = FALSE,
                     normalise = TRUE,
                     detailed_output = FALSE) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  ## Clean up the taxonomy
  data <- amp_rename(
    data = data,
    tax_class = tax_class,
    tax_empty = tax_empty,
    tax_level = tax_aggregate
  )

  # normalise counts
  if (isTRUE(normalise)) {
    data <- normaliseTo100(data)
  }

  # Aggregate to a specific taxonomic level
  abund1 <- aggregate_abund(
    abund = data$abund,
    tax = data$tax,
    tax_aggregate = tax_aggregate,
    tax_add = tax_add,
    calcSums = TRUE,
    format = "long"
  ) %>% 
    as.data.frame()

  ## Add group information
  suppressWarnings(
    if (group_by != "Sample") {
      if (length(group_by) > 1) {
        grp <- data.frame(Sample = rownames(data$metadata), Group = apply(data$metadata[, group_by], 1, paste, collapse = " "))
      } else {
        grp <- data.frame(Sample = rownames(data$metadata), Group = data$metadata[, group_by])
      }
      abund1$Group <- grp$Group[match(abund1$Sample, grp$Sample)]
      abund2 <- abund1
    } else {
      abund2 <- data.frame(abund1, Group = abund1$Sample)
    }
  )

  ## Take the average to group level
  abund3 <- data.table(abund2)[, Abundance := mean(Sum), by = list(Display, Group)] %>%
    setkey(Display, Group) %>%
    unique() %>%
    filter(Abundance > 0) %>%
    mutate(freq = 1)


  abund3$HA <- ifelse(abund3$Abundance > abund_thrh, 1, 0)
  temp3 <- group_by(abund3, Display) %>%
    summarise(Frequency = sum(freq), freq_A = sum(HA), Abundance = round(mean(Abundance), 2)) %>%
    as.data.frame()

  p <- ggplot(data = temp3, aes(x = Frequency, y = freq_A)) +
    xlab(paste("Observed in N ", group_by, "s", sep = "")) +
    theme_classic() +
    theme(
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "grey90")
    )


  if ((is.null(attributes(data)$normalised) | isFALSE(attributes(data)$normalised)) &
    isFALSE(normalise)) {
    p <- p + ylab(paste("Abundant in N ", group_by, "s (>", abund_thrh, "read counts)", sep = ""))
  } else {
    p <- p + ylab(paste("Abundant in N ", group_by, "s (>", abund_thrh, "%)", sep = ""))
  }

  if (plotly == F) {
    p <- p + geom_jitter(size = 3, alpha = 0.5)
  }

  if (isTRUE(plotly)) {
    data_plotly <- paste("Kingdom: ", data$tax[, 1], "<br>",
      "Phylum: ", data$tax[, 2], "<br>",
      "Class: ", data$tax[, 3], "<br>",
      "Order: ", data$tax[, 4], "<br>",
      "Family: ", data$tax[, 5], "<br>",
      "Genus: ", data$tax[, 6], "<br>",
      "Species: ", data$tax[, 7], "<br>",
      "OTU: ", data$tax[, 8],
      sep = ""
    )
    p <- p + geom_jitter(size = 2, alpha = 0.5, aes(text = data_plotly))
  }

  if (tax_aggregate == "OTU") {
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = data$tax, by = "OTU")
    temp3 <- core
  }

  if (isTRUE(plotly)) {
    ggplotly(p, tooltip = "text") %>%
      plotly::layout(showlegend = FALSE)
  } else if (isFALSE(plotly)) {
    if (detailed_output) {
      return(list(data = temp3, plot = p, abund = data$abund, tax = data$tax, metadata = data$metadata))
    } else if (!detailed_output) {
      return(p)
    }
  }
}
