#' Boxplot
#'
#' Generates boxplots of the most abundant taxa.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param order_group A vector to order the groups by.
#' @param order_y A vector to order the y-axis by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Genus"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{20})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_flip (\emph{logical}) Flip the axes of the plot axis. (\emph{default:} \code{FALSE})
#' @param plot_log (\emph{logical}) Log10-scale the plot. (\emph{default:} \code{FALSE})
#' @param adjust_zero Keep abundances of 0 in the calculation of medians by adding this value. (\emph{default:} \code{NULL})
#' @param point_size The size of points. (\emph{default:} \code{1})
#' @param sort_by Generic function name to use for sorting most abundant taxa, fx \code{mean}, \code{median}, or \code{sum}. (\emph{default:} \code{median})
#' @param plot_type Plot type. \code{"boxplot"} or \code{"point"}. (\emph{default:} \code{"boxplot"})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @import ggplot2
#' @importFrom dplyr arrange desc group_by summarise
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table setkey
#' @importFrom stats median
#'
#' @export
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # 10 boxplots grouped by WWTP with phylum name added
#' amp_boxplot(AalborgWWTPs,
#'   group_by = "Plant",
#'   tax_show = 10,
#'   tax_add = "Phylum"
#' )
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_boxplot <- function(data,
                        group_by = NULL,
                        sort_by = median,
                        plot_type = "boxplot",
                        point_size = 1,
                        tax_aggregate = "Genus",
                        tax_add = NULL,
                        tax_show = 20,
                        tax_empty = "best",
                        tax_class = NULL,
                        order_group = NULL,
                        order_y = NULL,
                        plot_flip = FALSE,
                        plot_log = FALSE,
                        adjust_zero = NULL,
                        normalise = TRUE) {
  
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
  
  # Group by sample if group_by is NULL, always coerce to factor
  sampleIDvarname <- colnames(data$metadata)[1] # also used later
  if(is.null(group_by)) {
    group_by <- sampleIDvarname
  }
  data$metadata[group_by] <- lapply(data$metadata[group_by], factor)
  
  # Aggregate to a specific taxonomic level and merge with chosen metadata group_by var(s)
  abund5 <- aggregate_abund(
    abund = data$abund,
    tax = data$tax,
    tax_aggregate = tax_aggregate,
    tax_add = tax_add,
    calcSums = TRUE,
    format = "long"
  ) %>%
    as.data.frame() %>% 
    merge(
      data.frame(
        Sample = data$metadata[[1]],
        .Group = apply(
          data$metadata[, group_by, drop = FALSE],
          1,
          paste,
          collapse = " "
        )
      ),
      by = "Sample"
    )
  
  ## Sort by chosen measure (median/mean/sum etc)
  TotalCounts <- abund5 %>% 
    group_by(Display) %>%
    summarise(measure = match.fun(sort_by)(Sum)) %>% 
    arrange(desc(measure))
  abund5$Display <- factor(abund5$Display, levels = rev(TotalCounts$Display))
  
  ## Subset to X most abundant levels
  if (is.numeric(tax_show)) {
    if (tax_show > nrow(TotalCounts)) {
      warning(paste0("There are only ", nrow(TotalCounts), " taxa, showing all"), call. = FALSE)
      tax_show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund5, Display %in% unique(TotalCounts$Display)[1:tax_show])
  } else if (!is.numeric(tax_show)) {
    tax_show <- as.character(tax_show)
    if (all(tolower(tax_show) == "all")) {
      abund7 <- abund5
    } else {
      abund7 <- filter(abund5, Display %in% tax_show)
    }
  }
  
  # filter returns a tibble in older versions
  abund7 <- as.data.frame(abund7)

  ## Add a small constant to handle ggplot2 removal of 0 values in log scaled plots
  if (!is.null(adjust_zero)) {
    abund7$Abundance[abund7$Abundance == 0] <- adjust_zero
  }

  ## Order y based on a vector
  if (length(order_y) > 1) {
    abund7$Display <- factor(abund7$Display, levels = order_y)
    abund7 <- subset(abund7, !is.na(Display))
  }

  ## plot the data
  if (group_by == sampleIDvarname) {
    p <- ggplot(abund7, aes(x = Display, y = Abundance))
  }
  if (group_by != sampleIDvarname) {
    if (!is.null(order_group)) {
      abund7$.Group <- factor(abund7$.Group, levels = rev(order_group))
    }
    p <- ggplot(abund7, aes(x = Display, y = Abundance, color = .Group))
  }

  p <- p +
    guides(col = guide_legend(reverse = TRUE)) +
    xlab("") +
    theme_classic() +
    theme(
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "grey90")
    )

  if ((is.null(attributes(data)$normalised) | isFALSE(attributes(data)$normalised)) &
    isFALSE(normalise)) {
    p <- p + ylab("Read counts")
  } else {
    p <- p + ylab("Relative Abundance (%)")
  }

  if (plot_flip == F) {
    p <- p + coord_flip()
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  if (plot_type == "point") {
    p <- p + geom_point(size = point_size)
  }
  if (plot_type == "boxplot") {
    p <- p + geom_boxplot(outlier.size = point_size)
  }
  if (plot_log == TRUE) {
    p <- p + scale_y_log10()
  }

  return(p)
}
