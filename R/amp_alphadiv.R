#' Alpha-diversity analysis
#'
#' Calculate alpha-diversity indices for each sample and combines with the metadata.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param measure Alpha-diversity measure(s) to be included if not all. A vector of one or more of:
#' \itemize{
#'   \item \code{"uniqueotus"}
#'   \item \code{"shannon"}
#'   \item \code{"simpson"}
#'   \item \code{"invsimpson"}
#' }
#' @param richness (\emph{logical}) Also calculate sample richness estimates (Chao1 and ACE) as calculated by \code{\link[vegan]{estimateR}}. (\emph{default:} \code{FALSE})
#' @param rarefy Rarefy species richness to this value before calculating alpha diversity and/or richness. Passed directly as the \code{sample} argument to \code{\link[vegan]{rrarefy}}. (\emph{default:} \code{NULL})
#' @param plot (\emph{logical}) Produce a boxplot instead of a table. (\emph{default:} \code{FALSE})
#' @param plot_scatter (\emph{logical}) If \code{TRUE} produce a scatter plot instead of a boxplot. (\emph{default:} \code{FALSE})
#' @param plot_group_by Group the samples by a categorical variable in the metadata. If \code{NULL} then all samples are shown in a single group.

#' @export
#' @importFrom dplyr arrange select
#' @importFrom vegan diversity estimateR
#' @importFrom data.table setDT melt
#'
#' @details The alpha-diversity indices are calculated per sample using the vegan function \code{\link[vegan]{diversity}}, where the read abundances are first rarefied using \code{\link[vegan]{rrarefy}} by the size of the \code{rarefy} argument. Refer to the vegan documentation for details about the different indices and how they are calculated. If no measure(s) are chosen, all diversity indices will be returned.
#'
#' @references
#' McMurdie, P.J. & Holmes, S. (2014). Waste not, want not: Why
#' rarefying microbiome data is inadmissible. \emph{PLoS Comput Biol}
#' \strong{10(4):} e1003531. DOI:\code{10.1371/journal.pcbi.1003531}
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Subsample/rarefy to 20000 reads and then calculate
#' # Shannon and Simpson alpha-diversity indices
#' alphadiversityresult <- amp_alphadiv(
#'   AalborgWWTPs,
#'   measure = c("shannon", "simpson"),
#'   rarefy = 20000
#' )
#'
#' # Explore the results in the data frame
#' # View(alphadiversityresult)
#' 
#' # Generate a plot instead
#' amp_alphadiv(
#'   AalborgWWTPs,
#'   measure = c("shannon", "simpson"),
#'   rarefy = 20000,
#'   plot = TRUE,
#'   plot_group_by = "Plant"
#' )
#' @return A data frame or a ggplot if \code{plot} is set to TRUE
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_alphadiv <- function(data,
                         measure = NULL,
                         richness = FALSE,
                         rarefy = NULL,
                         plot = FALSE,
                         plot_scatter = FALSE,
                         plot_group_by = NULL) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # what's the name of the first column (sample IDs)?
  sampleIDvarname <- colnames(data$metadata)[1]

  # check measures
  validMeasures <- c("uniqueotus", "shannon", "simpson", "invsimpson")
  if (is.null(measure)) {
    measure <- validMeasures
  } else if (!is.null(measure) & any(!measure %in% validMeasures)) {
    measure <- measure %>% tolower()
    warning("Some or none of the provided measures were not recognised, calculating all. Valid options are:\n", paste0(validMeasures, collapse = ", "), call. = FALSE)
    measure <- validMeasures
  }

  results <- data$metadata
  names <- results[[1]] # for making sure the ordering of the values in the vectors calculated later match the order of the metadata samples

  # nreads before rarefying
  RawReads <- colSums(data$abund)
  RawReads <- RawReads[names]
  results$RawReads <- RawReads

  if (!is.null(rarefy)) {
    data <- amp_rarefy(data, rarefy)
  }

  # Add Reads column
  Reads <- colSums(data$abund)
  Reads <- Reads[names]
  results$Reads <- Reads

  # warning from phyloseq::estimate_richness
  if (!any(data$abund == 1)) {
    warning("The data you have provided does not have\n",
      "any singletons. This is highly suspicious. Results of richness\n",
      "estimates (for example) are probably unreliable, or wrong, if you have already\n",
      "trimmed low-abundance taxa from the data.\n", "\n",
      "We recommend that you find the untrimmed data and retry.",
      call. = FALSE
    )
  }

  tabund <- t(data$abund)
  if (any("uniqueotus" %in% measure) | is.null(measure)) {
    uniqueOTUs <- colSums(data$abund > 0) # not transposed
    results$uniqueOTUs <- uniqueOTUs[names]
  }
  if (any("shannon" %in% measure)) {
    Shannon <- vegan::diversity(tabund, index = "shannon")
    # For some reason vegan throws away names if only 1 sample
    if (length(Shannon) == 1) {
      results$Shannon <- Shannon
    } else if (length(Shannon) > 1) {
      results$Shannon <- Shannon[names]
    }
  }
  if (any("simpson" %in% measure)) {
    Simpson <- vegan::diversity(tabund, index = "simpson")
    # For some reason vegan throws away names if only 1 sample
    if (length(Simpson) == 1) {
      results$Simpson <- Simpson
    } else if (length(Simpson) > 1) {
      results$Simpson <- Simpson[names]
    }
  }
  if (any("invsimpson" %in% measure)) {
    invSimpson <- vegan::diversity(tabund, index = "invsimpson")
    # For some reason vegan throws away names if only 1 sample
    if (length(invSimpson) == 1) {
      results$invSimpson <- invSimpson
    } else if (length(invSimpson) > 1) {
      results$invSimpson <- invSimpson[names]
    }
  }
  if (isTRUE(richness)) {
    richness <- t(vegan::estimateR(tabund)) %>% as.data.frame()
    richness <- richness[names, , drop = FALSE]
    results$Chao1 <- richness[, "S.chao1"]
    results$ACE <- richness[, "S.ACE"]
  }
  
  if (isTRUE(plot)) {
    setDT(results)
    plot_data <- melt(
      results,
      id.vars = c(sampleIDvarname, plot_group_by),
      measure.vars = colnames(results)[tolower(colnames(results)) %in% tolower(measure)],
      variable.name = "measure",
      value.name = "value"
    )
    ggplot(
      plot_data,
      aes(
        x = if(is.null(plot_group_by)) {
          sampleIDvarname
        } else {
          eval(parse(text = plot_group_by))
        },
        y = value
      )
    ) +
    {
      if (!isTRUE(plot_scatter)) {
        geom_boxplot()
      } else if (isTRUE(plot_scatter)) {
        geom_point()
      }
    } +
    scale_fill_viridis_d(
      option = "turbo"
    ) +
    facet_wrap(~measure, nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(
      strip.background = element_rect(
        colour = "gray8", fill = "gray98"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.major.y = element_line(color = "gray70"),
      strip.text = element_text(size = 16),
      axis.ticks.x = element_blank(),
      axis.text.x = if(is.null(plot_group_by)) {
        element_blank()
      } else {
        element_text(angle = 45, size = 14, hjust = 1, vjust = 1)
      },
      axis.text.y = element_text(size = 14),
      axis.title = element_blank()
    )

  } else {
    # arrange by RawReads no matter if rarefied or not, remove the column if rarefy is not set
    results <- results %>%
      dplyr::arrange(RawReads) %>%
      {
        if (is.null(rarefy)) dplyr::select(., -RawReads) else .
      }
    return(results)
  }
}

#' @rdname amp_alphadiv
#' @export
amp_alpha_diversity <- amp_alphadiv
