#' Rarefaction curve
#'
#' Generates a rarefaction curve (number of reads vs number of observed OTUs) for each sample.
#'
#' @usage amp_rarecurve(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param stepsize Step size for the curves. Lower is prettier but takes more time to generate. (\emph{default:} \code{1000})
#' @param color_by Color curves by a variable in the metadata. (\emph{default:} \code{NULL})
#' @param facet_by Split the plot into subplots based on a variable in the metadata. (\emph{default:} \code{NULL})
#' @param facet_scales If \code{facet_by} is set, should the axis scales of each subplot be fixed (\code{fixed}), free (\code{"free"}), or free in one dimension (\code{"free_x"} or \code{"free_y"})? (\emph{default:} \code{"fixed"})
#'
#' @export
#'
#' @import ggplot2
#' @importFrom vegan rarefy
#' @importFrom plyr ldply
#'
#' @return A ggplot2 object.
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Rarecurve
#' amp_rarecurve(AalborgWWTPs, facet_by = "Plant")
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_rarecurve <- function(data,
                          stepsize = 1000,
                          color_by = NULL,
                          facet_by = NULL,
                          facet_scales = "fixed") {
  if (class(data) != "ampvis2") {
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  }

  maxreads <- max(colSums(data$abund))
  if (maxreads < stepsize) {
    stop("\"stepsize\" too high, maximum number of reads in any sample is: ", maxreads, call. = FALSE)
  }

  abund <- data[["abund"]] %>%
    as.matrix() %>%
    t()
  if (!identical(all.equal(abund, round(abund)), TRUE)) {
    stop("Function accepts only integers (counts)", call. = FALSE)
  }

  tot <- rowSums(abund)
  nr <- nrow(abund)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = stepsize)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    drop(vegan::rarefy(abund[i, ], n))
  })
  names(out) <- names(tot)
  df <- plyr::ldply(out, function(x) {
    data.frame(
      Species = x,
      Reads = attr(x, "Subsample"),
      check.names = FALSE
    )
  },
  .id = colnames(data[["metadata"]])[1]
  )

  gg <- merge(data[["metadata"]],
    df,
    by = 1
  )
  p <- ggplot(
    gg,
    aes_string(
      x = "Reads",
      y = "Species",
      group = colnames(data[["metadata"]])[1],
      color = color_by
    )
  ) +
    geom_line() +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, vjust = 0.3, angle = 90),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    xlab("Sequencing depth (reads)") +
    ylab("Number of observed OTUs")

  if (!is.null(facet_by)) {
    p <- p +
      facet_wrap(reformulate(facet_by),
        scales = facet_scales
      ) +
      theme(
        strip.background = element_rect(colour = NA, fill = "grey95"),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position = "bottom",
        strip.text = element_text(size = 10),
        legend.title = element_blank()
      )
  }
  return(p)
}
