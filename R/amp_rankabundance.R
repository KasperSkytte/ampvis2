#' @title Rank abundance plot
#'
#' @description Generates a rank abundance curve (rank abundance vs cumulative read abundance), optionally with standard deviation from mean intervals.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by (\emph{required}) Group the samples by a variable in the metadata.
#' @param showSD Show standard deviation from mean or not. (\emph{default:} \code{TRUE})
#' @param log10_x Log10-transform the x-axis or not. Often most variation is observed among the most abundant OTU's, log10-transforming the x-axis will highlight this better. (\emph{default:} \code{TRUE})
#'
#' @return A ggplot2 object.
#' @import ggplot2
#' @importFrom stats sd
#'
#' @export
#'
#' @details Currently only OTU level is supported.
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Rank abundance plot
#' amp_rankabundance(AalborgWWTPs, group_by = "Plant", showSD = TRUE, log10_x = TRUE)
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_rankabundance <- function(data,
                              group_by,
                              showSD = TRUE,
                              log10_x = TRUE) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # melt
  d <- amp_export_long(
    normaliseTo100(data),
    metadata_vars = group_by,
    tax_levels = "OTU"
  )

  # group up and summarise
  d[, groupSum := sum(count), by = c("OTU", group_by)]
  setorderv(d, c(group_by, "groupSum"), order = -1)
  d[, prop := groupSum / sum(groupSum) * 100, by = group_by]
  d[, cumProp := cumsum(prop), by = group_by]
  d[, rank := as.integer(!duplicated(groupSum)), by = c("OTU", group_by)]
  d[, rank := cumsum(rank), by = group_by]
  d[, sd := sd(count, na.rm = TRUE), by = c("OTU", group_by)]
  d <- d[, .SD[.N], by = c("OTU", group_by)]

  # generate plot
  plot <- ggplot(
    d,
    aes_string(
      x = "rank",
      y = "cumProp",
      group = group_by,
      color = group_by
    )
  ) +
    geom_line(size = 1) +
    ylim(0, 100) +
    xlab("Rank abundance") +
    ylab("Cumulative read abundance (%)") +
    theme_classic() +
    {
      if (isTRUE(showSD)) {
        geom_ribbon(
          aes_string(
            ymin = "cumProp - sd",
            ymax = "cumProp + sd",
            fill = group_by
          ),
          alpha = 0.2,
          size = 0
        )
      }
    } +
    {
      if (isTRUE(log10_x)) {
        scale_x_log10()
      }
    }
  return(plot)
}

#' @rdname amp_rankabundance
#' @export
amp_rank_abundance <- amp_rankabundance
