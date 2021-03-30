#' @title Grouped core community analysis
#'
#' @description Core community plot to investigate how many of the most abundant OTU's comprise the "core" of groups of samples, their abundances etc.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by (\emph{required}) Character vector of one or more variable name(s) in the sample metadata which contain the desired grouping of samples, fx where the samples have been taken. Can also be done by sample, just provide the name of the first variable containing unique sample ID's if so.
#' @param core_pct Threshold in percent for defining "abundant" a.k.a. "core" OTU's for each group. (\emph{default:} \code{80})
#' @param margin_plots Character vector defining which margin plots to show. Margin plots display the cumulative read abundances of all OTU's per group sharing either x or y axis with the main plot. Any of:
#' \itemize{
#'    \item \code{"x"}: Only show x axis margin plot.
#'    \item \code{"y"}: Only show y axis margin plot.
#'    \item \code{"xy"} (\emph{default}): Show both x and y margin plots.
#'    \item \code{""} or \code{NULL} (or anything else): Don't show any margin plots.
#'    }
#' @param margin_plot_values_size The size of the values indicated in the margin plots on top of the bars. Set to \code{0} to disable. (\emph{default:} \code{3})
#' @param widths Numeric vector with relative widths of the main and y margin plots. (\emph{default:} \code{c(5,1)})
#' @param heights Numeric vector with relative widths of the main and x margin plots. (\emph{default:} \code{c(1,5)})
#'
#' @details This analysis only makes sense without aggregating OTU's to any taxonomic level, or else it will be biased by taxonomy and only be done on OTU's that have been classified, which is rarely all.
#' @section Saving plot with \code{\link[ggplot2]{ggsave}}:
#' When any margin plots are generated \code{\link{amp_core}} returns a list of ggplot objects to allow adjusting themes etc. of the individual subplots. The list is of class \code{coreplot} and a matching print function for the S3 class then stitches together the individual plots using the \code{\link{patchwork}} package. Therefore to save the plot with \code{\link[ggplot2]{ggsave}} simply pass on the plot object explicitly and wrap it in print(), see examples. This is not necessary if no margin plots are generated, as the returned object is then a regular ggplot object.
#'
#' @return If no margin plots a \code{ggplot} object, otherwise a list with ggplot objects.
#'
#' @import ggplot2
#' @importFrom data.table uniqueN setorderv fifelse
#'
#' @export
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @references
#' Marta Nierychlo, Kasper Skytte Andersen, Yijuan Xu, Nicholas Green, Chenjing Jiang, Mads Albertsen, Morten Simonsen Dueholm, Per Halkjær Nielsen (2020): "MiDAS 3: An ecosystem-specific reference database, taxonomy and knowledge platform for activated sludge and anaerobic digesters reveals species-level microbiome composition of activated sludge", Water Research, Volume 182. \url{https://doi.org/10.1016/j.watres.2020.115955}
#' Saunders, Aaron M; Albertsen, Mads; Vollertsen, Jes; Nielsen, Per H (2016): "The activated sludge ecosystem contains a core community of abundant organisms", ISME Journal 10, 11-20. \url{https://doi.org/10.1038/ismej.2015.117}
#'
#' @examples
#' # load example data
#' data("MiDAS")
#'
#' # generate core analysis plot grouped by WWTP
#' plot <- amp_core(
#'   data = MiDAS,
#'   group_by = "Plant",
#'   core_pct = 80,
#'   margin_plots = "xy",
#'   margin_plot_values_size = 0 # set to 0 to not show values in margin plots
#' )
#'
#' # adjust axes manually
#' plot$mainplot <- plot$mainplot +
#'   scale_x_discrete(breaks = seq(from = 0, to = 36, by = 5)) +
#'   scale_y_discrete(breaks = seq(from = 0, to = 36, by = 5))
#' plot$marginplot_x <- plot$marginplot_x + scale_x_discrete(breaks = seq(from = 0, to = 36, by = 5))
#' plot$marginplot_y <- plot$marginplot_y + scale_x_discrete(breaks = seq(from = 0, to = 36, by = 5))
#'
#' # show plot
#' plot
#'
#' # To save the plot with ggsave() wrap the plot object in print()
#' # ggsave("plot.png", print(plot))
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_core <- function(data,
                     group_by,
                     core_pct = 80,
                     margin_plots = "xy",
                     margin_plot_values_size = 3,
                     widths = c(5, 1),
                     heights = c(1, 5)) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # melt into long format
  d <- amp_export_long(
    data,
    metadata_vars = group_by,
    tax_levels = "OTU"
  )

  # grouped calculations and ordering
  gg <- d[, .(sum = sum(count)), by = c("OTU", group_by)]
  setorderv(gg, c(group_by, "sum"), order = -1)
  gg[, OTUprop := sum / sum(sum) * 100, by = group_by]
  gg[, cumOTUprop := cumsum(OTUprop), by = group_by]
  gg[, core := fifelse(cumOTUprop < ..core_pct, TRUE, FALSE)]
  gg[, nObs := sum(sum > 0), by = OTU]
  gg[, nCore := sum(core), by = OTU]
  gg[, totalReads := sum(sum)]

  # create list of plots to return and print using a custom print method
  # (to allow adjusting plots individually)
  outlist <- structure(
    .Data = list(
      mainplot = NULL,
      marginplot_x = NULL,
      marginplot_y = NULL,
      data = gg
    ),
    class = "coreplot",
    widths = widths,
    heights = heights
  )

  # generate main plot
  outlist$mainplot <- ggplot(
    # summarise data (reducing size)
    gg[
      ,
      .(nOTUs = uniqueN(OTU)),
      by = .(nObs, nCore)
    ],
    aes(
      x = as.factor(nObs), # factors to align correctly with margin plots
      y = as.factor(nCore), # factors to align correctly with margin plots
      color = nOTUs,
      size = nOTUs
    )
  ) +
    geom_point() +
    xlab(paste0("Observed in N groups/samples (", group_by, ")")) +
    ylab(paste0("Part of top ", core_pct, "% of all reads\n in N groups/samples")) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(5, "RdYlBu")), trans = "log10") +
    theme_minimal()

  # generate x margin plot if chosen
  if (any(tolower(margin_plots) %in% c("x", "xy", "yx"))) {
    outlist$marginplot_x <- ggplot(
      gg[
        ,
        .(nObsSum = sum(sum) / unique(totalReads) * 100),
        by = .(nObs)
      ],
      aes(as.factor(nObs), nObsSum)
    ) +
      geom_col() +
      {
        if (margin_plot_values_size != 0) {
          geom_text(
            aes(
              label = round(nObsSum, 1)
            ),
            vjust = -1,
            size = margin_plot_values_size
          )
        }
      } +
      ylab("Cumulative OTU \nread abundance (%)") +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
  }

  # generate y margin plot if chosen
  if (any(tolower(margin_plots) %in% c("y", "xy", "yx"))) {
    outlist$marginplot_y <- ggplot(
      gg[
        ,
        .(nCoreSum = sum(sum) / unique(totalReads) * 100),
        by = .(nCore)
      ],
      aes(
        x = as.factor(nCore),
        y = nCoreSum
      )
    ) +
      geom_col() +
      {
        if (margin_plot_values_size != 0) {
          geom_text(
            aes(
              label = round(nCoreSum, 1)
            ),
            vjust = -1,
            size = margin_plot_values_size
          )
        }
      } +
      ylab("Cumulative OTU \nread abundance (%)") +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      coord_flip()
  }

  # return main plot only if no margin plots
  if (is.null(outlist$marginplot_x) && is.null(outlist$marginplot_y)) {
    return(outlist$mainplot)
  } else {
    return(outlist)
  }
}
