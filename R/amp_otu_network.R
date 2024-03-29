#' Network plot
#'
#' Generates network plot of taxa and samples based on ggnet2.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param color_by A metadata variable to color the samples by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Phylum"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{10})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param min_abundance Minimum taxa abundance pr. sample. (\emph{default:} \code{0})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param ... Additional arguments passed on to \code{\link[GGally]{ggnet2}}.
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @import ggplot2
#' @importFrom dplyr filter arrange group_by mutate select summarise desc
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table data.table setkey
#'
#' @export
#'
#' @section Preserving relative abundances in a subset of larger data:
#' See \code{?\link{amp_filter_samples}} or the \href{https://kasperskytte.github.io/ampvis2/articles/faq.html#preserving-relative-abundances-in-a-subset-of-larger-data}{ampvis2 FAQ}.
#'
#' @details See \code{\link[GGally]{ggnet2}}
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # OTU network plot
#' amp_otu_network(AalborgWWTPs)
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_otu_network <- function(data,
                            min_abundance = 0,
                            color_by = NULL,
                            tax_aggregate = "Phylum",
                            tax_add = NULL,
                            tax_show = 10,
                            tax_class = NULL,
                            tax_empty = "best",
                            normalise = TRUE,
                            ...) {
  checkReqPkg("GGally")
  checkReqPkg("network")
  checkReqPkg("sna")

  ### Data must be in ampvis2 format
  is_ampvis2(data)

  ## Clean up the taxonomy
  data <- amp_rename(
    data = data,
    tax_class = tax_class,
    tax_empty = tax_empty,
    tax_level = tax_aggregate
  )

  ## SampleID column is used to merge data later, so it must be there!
  colnames(data$metadata)[1] <- "SampleID"

  # normalise counts
  if (isTRUE(normalise)) {
    data <- normaliseTo100(data)
  }

  # Aggregate to a specific taxonomic level
  abund3 <- aggregate_abund(
    abund = data$abund,
    tax = data$tax,
    tax_aggregate = tax_aggregate,
    tax_add = tax_add,
    calcSums = TRUE,
    format = "long"
  ) %>%
    as.data.frame()

  ## Add group information
  abund5 <- data.frame(abund3, .Group = abund3$Sample)

  ## Take the average to group level

  abund6 <- data.table(abund5)[, Abundance := mean(Sum), by = list(Display, .Group)] %>%
    setkey(Display, .Group) %>%
    unique() %>%
    as.data.frame() %>%
    select(-Sum)

  ## Find the X most abundant levels
  TotalCounts <- group_by(abund6, Display) %>%
    summarise(Abundance = sum(Abundance)) %>%
    arrange(desc(Abundance))

  if (tax_show > nrow(TotalCounts)) {
    tax_show <- nrow(TotalCounts)
  }
  abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax_show])

  ## Convert to network
  netw <- data.frame(
    SampleID = as.character(abund7$.Group),
    Taxa = abund7$Display,
    Abundance = abund7$Abundance,
    stringsAsFactors = F
  ) %>%
    subset(Abundance > min_abundance) %>%
    # Subset to dominant species in each sample
    select(-Abundance) %>%
    network(directed = FALSE)

  ## Add data to nodes
  x <- data.frame(SampleID = network.vertex.names(netw), stringsAsFactors = F)

  xsamples <- filter(x, !grepl("Taxa", SampleID)) %>%
    merge(data$metadata, all.x = T, by = 1)

  if (is.null(color_by)) {
    xsamples$Description <- "Sample"
  } else {
    xsamples$Description <- as.character(xsamples[, color_by])
  }

  set.vertex.attribute(netw, "snames", c(rep("Sample", nrow(xsamples)), rep("Taxa", nrow(x) - nrow(xsamples))))
  set.vertex.attribute(netw, "stype", c(xsamples$Description, rep("Taxa", nrow(x) - nrow(xsamples))))
  set.vertex.attribute(netw, "nsize", c(rep(3, nrow(xsamples)), rep(1, nrow(x) - nrow(xsamples))))

  ## Make network plot

  p <- ggnet2(netw,
    color = "stype",
    node.alpha = 0.7,
    edge.color = "grey80",
    node.size = "nsize",
    ...
  ) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_size_discrete(guide = F, range = c(3, 6))

  return(p)
}
