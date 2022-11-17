#' Heatmap
#'
#' Generates a heatmap of amplicon data by using sample metadata to aggregate samples and taxonomy to aggregate OTUs.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by (\emph{recommended}) Group the samples by a categorical variable in the metadata. If \code{NULL} then all samples are shown.
#' @param facet_by Facet the samples by a categorical variable in the metadata.
#' @param scale_by Scale the abundances by a variable in the metadata.
#' @param normalise_by A variable or a specific sample in the metadata to normalise the counts by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Phylum"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{10})
#' @param showRemainingTaxa Add an additional row at the bottom displaying the sum of all remaining taxa that are not part of the top \code{tax_show} most abundant taxa. (\emph{default:} \code{FALSE})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param measure Calculate and display either \code{"mean"}, \code{"max"} or \code{"median"} across the groups. (\emph{default:} \code{"mean"})
#' @param sort_by Sorts the heatmap by the most abundant taxa in a specific sample or group of samples. Provide a sample name or a specific value of the group defined by the \code{"group_by"} argument, e.g. \code{"Treatment A"}.
#' @param order_x_by Reorder the x axis by providing a character vector of the x axis values in the desired order, or \code{"cluster"} for hierarchical clustering by \code{\link[stats]{hclust}}.
#' @param order_y_by Reorder the y axis by providing a character vector of the y axis values in the desired order, or \code{"cluster"} for hierarchical clustering by \code{\link[stats]{hclust}}.
#' @param plot_values (\emph{logical}) Plot the values on the heatmap or not. (\emph{default:} \code{TRUE})
#' @param plot_values_size The size of the plotted values. (\emph{default:} \code{4})
#' @param plot_legendbreaks A vector of breaks for the abundance legend, fx \code{c(1, 10, 20)}.
#' @param plot_colorscale The type of scale used for the coloring of abundances, either \code{"sqrt"} or \code{"log10"}. (\emph{default:} \code{"log10"})
#' @param plot_na (\emph{logical}) Whether to color missing values with the lowest color in the scale or not. (\emph{default:} \code{TRUE})
#' @param min_abundance All values below this value are given the same color. (\emph{default:} \code{0.1})
#' @param max_abundance All values above this value are given the same color.
#' @param color_vector Vector of colors for the colorscale, e.g. \code{c("white", "red")}.
#' @param round Number of digits to show with the values. (\emph{default:} \code{1})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param textmap (\emph{logical}) Return a data frame to print as raw text instead of a ggplot2 object. (\emph{default:} \code{FALSE})
#' @param plot_functions Return a 2-column grid plot instead, showing known functional information about the Genus-level OTUs next to the heatmap. By default, this functional information is retrieved directly from \url{midasfieldguide.org}. When using this feature, make sure that either \code{tax_aggregate} or \code{tax_add} is set to "Genus" and that Genus is the lowest level in either. (\emph{default:} \code{FALSE})
#' @param function_data If \code{plot_functions} is set to \code{TRUE}: A data frame with functional information at Genus level. The first column must be the Genus names and any other column(s) can be any property or metabolic function of the individual Genera. If \code{NULL} then data will be retrieved directly from \url{midasfieldguide.org}. (\emph{default:} \code{NULL})
#' @param functions If \code{plot_functions} is set to \code{TRUE}: A vector with the functions to be displayed (column names in the \code{functions_data} data frame). If data is succesfully retrieved from \url{midasfieldguide.org} then available functions can be listed with \code{colnames(.ampvis2_midasfg_function_data)}. (\emph{default:} \code{c("MiDAS","Filamentous", "AOB", "NOB", "PAO", "GAO")})
#' @param rel_widths If \code{plot_functions} is set to \code{TRUE}: A vector with the relative widths of the heatmap and function grid when \code{plot_functions = TRUE}. (\emph{default:} \code{c(0.75, 0.25)})
#'
#' @return A ggplot2 object, or a list of ggplot2 objects if \code{plot_functions = TRUE}. A data frame if \code{textmap = TRUE}.
#'
#' @export
#'
#' @section Preserving relative abundances in a subset of larger data:
#' By default the raw read counts in the abundance matrix are normalised (transformed to percentages) by some plotting functions automatically (for example \code{\link{amp_heatmap}}, \code{\link{amp_timeseries}}, and more). This means that the relative abundances shown will be calculated based on the remaining taxa after the subset, not including the removed taxa, if any. To circumvent this, set \code{normalise = TRUE} when subsetting with the \code{\link{amp_filter_taxa}} and \code{\link{amp_filter_samples}} functions, and then set \code{normalise = FALSE} in the plotting function. This will transform the OTU counts to relative abundances BEFORE the subset, and setting \code{normalise = FALSE} will skip the transformation in the plotting function, see the example below.
#'
#' \preformatted{
#' data("MiDAS")
#' subsettedData <- amp_filter_samples(MiDAS,
#'                                     Plant \%in\% c("Aalborg West", "Aalborg East"),
#'                                     normalise = TRUE
#'                                     )
#' amp_heatmap(subsettedData,
#'             group_by = "Plant",
#'             tax_aggregate = "Phylum",
#'             tax_add = "Genus",
#'             normalise = FALSE
#'             )
#' }
#'
#' @section Saving plot with \code{\link[ggplot2]{ggsave}}:
#' When \code{plot_functions = TRUE} a list of ggplot objects is returned to allow adjusting themes etc. of the individual subplots. The list is of class \code{hmfunplot} and a matching print function for the S3 class then stitches together the individual plots using the \href{https://patchwork.data-imaginist.com/}{patchwork} package. Therefore to save the plot with \code{\link[ggplot2]{ggsave}} simply pass on the plot object explicitly and wrap it in print(), see examples. This is not necessary if \code{plot_functions = FALSE}, as the returned object is then a regular ggplot object.
#'
#' @section Accessing detailed raw data:
#' The complete raw data used to generate any ggplot can always be accessed with \code{ggplot2_object$data} when the plot is saved as a ggplot2 object. Additionally, a "textmap" version of the generated heatmap can also be generated by setting \code{textmap = TRUE} to only extract the raw data as shown on the particular heatmap, see examples.
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Heatmap grouped by WWTP
#' amp_heatmap(AalborgWWTPs, group_by = "Plant")
#'
#' # Heatmap of 20 most abundant Genera (by mean) grouped by WWTP, split by Year,
#' # values not plotted for visibility, phylum name added, colorscale adjusted manually,
#' # and show the sum of remaining taxa not part of the top 20 most abundant taxa
#' amp_heatmap(AalborgWWTPs,
#'   group_by = "Plant",
#'   facet_by = "Year",
#'   plot_values = FALSE,
#'   tax_show = 20,
#'   showRemainingTaxa = TRUE,
#'   tax_aggregate = "Genus",
#'   tax_add = "Phylum",
#'   color_vector = c("white", "red"),
#'   plot_colorscale = "sqrt",
#'   plot_legendbreaks = c(1, 5, 10)
#' )
#'
#' # Heatmap with known functional traits about the Genera shown to the right
#' # By default this information is retrieved directly from midasfieldguide.org
#' # but you can provide your own with the function_data argument as shown with
#' # the textmap further down
#'
#' suppressWarnings(
#'   heatmapwfunctions <- amp_heatmap(AalborgWWTPs,
#'     group_by = "Plant",
#'     tax_aggregate = "Genus",
#'     plot_functions = TRUE,
#'     functions = c("PAO", "GAO", "AOB", "NOB")
#'   )
#' )
#'
#' class(heatmapwfunctions)
#'
#' # To save the plot with ggsave() wrap the plot object in print()
#' # ggsave("plot.png", print(heatmapwfunctions))
#' # The special class is essentially just a list of ggplots, allowing you to 
#' # customize them individually with standard ggplot functions.
#'
#' # A raw text version of the heatmap can be printed or saved as a data frame with textmap = TRUE.
#' textmap <- amp_heatmap(AalborgWWTPs,
#'   group_by = "Plant",
#'   tax_aggregate = "Genus",
#'   plot_functions = TRUE,
#'   function_data = midasfunctions,
#'   functions = c("PAO", "GAO", "AOB", "NOB"),
#'   textmap = TRUE
#' )
#' textmap
#' @import ggplot2
#' @importFrom dplyr filter desc arrange group_by mutate summarise
#' @importFrom tidyr gather spread
#' @importFrom data.table as.data.table data.table setkey dcast melt setDT setDF rbindlist
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales squish
#' @importFrom stats dist hclust median reformulate
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_heatmap <- function(data,
                        group_by = NULL,
                        facet_by = NULL,
                        normalise = TRUE,
                        tax_aggregate = "Phylum",
                        tax_add = NULL,
                        tax_show = 10,
                        showRemainingTaxa = FALSE,
                        tax_class = NULL,
                        tax_empty = "best",
                        order_x_by = NULL,
                        order_y_by = NULL,
                        plot_values = TRUE,
                        plot_values_size = 4,
                        plot_legendbreaks = NULL,
                        plot_colorscale = "log10",
                        plot_na = TRUE,
                        measure = "mean",
                        min_abundance = 0.1,
                        max_abundance = NULL,
                        sort_by = NULL,
                        normalise_by = NULL,
                        scale_by = NULL,
                        color_vector = NULL,
                        round = 1,
                        textmap = FALSE,
                        plot_functions = FALSE,
                        function_data = NULL,
                        functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
                        rel_widths = c(0.75, 0.25)) {

  ### Data must be in ampvis2 format
  is_ampvis2(data)

  ## Clean up the taxonomy
  data <- amp_rename(
    data = data,
    tax_class = tax_class,
    tax_empty = tax_empty,
    tax_level = tax_aggregate
  )

  # tax_add and tax_aggregate can't be the same
  if (!is.null(tax_aggregate) & !is.null(tax_add)) {
    if (identical(tax_aggregate, tax_add)) {
      stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
    }
  }

  # Checks an data if plot_functions = TRUE
  if (isTRUE(plot_functions)) {
    if (!any("Genus" %in% c(tax_add, tax_aggregate))) {
      stop("One of the arguments tax_add or tax_aggregate must contain \"Genus\"", call. = FALSE)
    }

    # Retrieve the function data from MiDAS fieldguide if none is provided
    # Only once per session, save in a hidden object .ampvis2_midasfg_function_data
    if (is.null(function_data)) {
      if (!exists(".ampvis2_midasfg_function_data", envir = .GlobalEnv)) {
        function_data <- try(getMiDASFGData(), silent = TRUE)
        if (inherits(function_data, "try-error")) {
          warning("Can't reach the midasfieldguide.org API to download functional data just now. The reason can be issues with either the site or your internet connection. Using the data set \"midasfunctions\" instead, which is probably not up-to-date. You can also supply your own data frame by using the function_data argument.", call. = FALSE)
          function_data <- midasfunctions
        } else {
          function_data <- extractFunctions(function_data)
        }
        assign(".ampvis2_midasfg_function_data", function_data, envir = .GlobalEnv)
      } else if (exists(".ampvis2_midasfg_function_data", envir = .GlobalEnv)) {
        function_data <- get(".ampvis2_midasfg_function_data", envir = .GlobalEnv)
      }
    }

    # check if the chosen functions are in the function_data
    checkFuncs <- functions %in% colnames(function_data)
    if (!all(checkFuncs)) {
      stop(
        paste0(
          "\"",
          paste0(functions[!checkFuncs], collapse = "\", \""),
          "\" not found in function_data. Available ones are:\n",
          paste0(colnames(function_data)[-1], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  ## Coerce the group_by and facet_by variables to factor to always be considered categorical.
  # Fx Year is automatically loaded as numeric by R, but it should be considered categorical.
  ## Grouping a heatmap by a continuous variable doesn't make sense
  if (!is.null(group_by)) {
    data$metadata[group_by] <- lapply(data$metadata[group_by], factor)
  }

  if (!is.null(facet_by)) {
    if (is.null(group_by)) {
      group_by <- names(data$metadata)[[1]]
    }
    data$metadata[facet_by] <- lapply(data$metadata[facet_by], factor)
  }

  ## Scale the data by a selected metadata sample variable
  if (!is.null(scale_by)) {
    variable <- as.numeric(data$metadata[, scale_by])
    data$abund <- t(t(data$abund) * variable)
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
  if (!is.null(facet_by)) {
    ogroup <- group_by
    group_by <- c(group_by, facet_by)
  }

  suppressWarnings(
    if (!is.null(group_by)) {
      if (length(group_by) > 1) {
        grp <- data.frame(Sample = data$metadata[, 1], .Group = apply(data$metadata[, group_by], 1, paste, collapse = " "))
        oldGroup <- unique(cbind.data.frame(data$metadata[, group_by], .Group = grp$.Group))
      } else {
        grp <- data.frame(Sample = data$metadata[, 1], .Group = data$metadata[, group_by])
      }
      abund3$.Group <- grp$.Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else {
      abund5 <- data.frame(abund3, .Group = abund3$Sample)
    }
  )

  ## Take the average to group level
  if (measure == "mean") {
    abund6 <- data.table(abund5)[, Abundance := mean(Sum), by = list(Display, .Group)] %>%
      setkey(Display, .Group) %>%
      unique() %>%
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = sum(Abundance)) %>%
      arrange(desc(Abundance))
  }

  if (measure == "max") {
    abund6 <- data.table(abund5)[, Abundance := max(Sum), by = list(Display, .Group)] %>%
      setkey(Display, .Group) %>%
      unique() %>%
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = max(Abundance)) %>%
      arrange(desc(Abundance))
  }

  if (measure == "median") {
    abund6 <- data.table(abund5)[, Abundance := median(Sum), by = list(Display, .Group)] %>%
      setkey(Display, .Group) %>%
      unique() %>%
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = median(Abundance)) %>%
      arrange(desc(Abundance))
  }

  if (!is.null(sort_by)) {
    if (is.null(group_by)) {
      if (!any(sort_by %in% abund6$Sample)) {
        stop("Can't find \"", sort_by, "\" among sample names", call. = FALSE)
      }
      TotalCounts <- filter(abund6, Sample == sort_by) %>%
        arrange(desc(Abundance))
    } else if (!is.null(group_by)) {
      if (any(sort_by %in% abund6$.Group) & any(sort_by %in% abund6$Sample)) {
        stop(paste0(sort_by, " is both found among samples and in the group_by variable (", group_by, "). Cannot sort by both a sample and a group."), call. = FALSE)
      } else if (!any(sort_by %in% abund6$.Group) & any(sort_by %in% abund6$Sample)) {
        TotalCounts <- filter(abund6, Sample == sort_by) %>%
          arrange(desc(Abundance))
      } else if (any(sort_by %in% abund6$.Group) & !any(sort_by %in% abund6$Sample)) {
        TotalCounts <- filter(abund6, .Group == sort_by) %>%
          arrange(desc(Abundance))
      } else if (!any(sort_by %in% abund6$.Group) & !any(sort_by %in% abund6$Sample)) {
        stop("Can't find \"", sort_by, "\" among sample or group names", call. = FALSE)
      }
    }
  }

  ## Subset to X most abundant levels
  if (is.numeric(tax_show)) {
    if (tax_show > nrow(TotalCounts)) {
      warning(paste0("There are only ", nrow(TotalCounts), " taxa, showing all"), call. = FALSE)
      tax_show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund6, Display %in% unique(TotalCounts$Display)[1:tax_show])
  } else if (!is.numeric(tax_show)) {
    tax_show <- as.character(tax_show)
    if (all(tolower(tax_show) == "all")) {
      abund7 <- abund6
    } else {
      abund7 <- filter(abund6, Display %in% tax_show)
    }
  }
  abund7 <- as.data.frame(abund7)

  ## Normalise to a specific group or sample (The Abundance of the group is set as 1)
  if (!is.null(normalise_by)) {
    if (!normalise_by %chin% unique(abund7$.Group)) {
      stop(paste0(normalise_by, " is not found among group names, cannot normalise"), call. = FALSE)
    }
    setDT(abund7)
    normalise_by <- abund7[.Group == normalise_by, Abundance]
    abund7[, Abundance := Abundance / normalise_by, by = .Group]
    setDF(abund7)
  }

  ## Order.y
  if (is.null(order_y_by)) {
    abund7$Display <- factor(abund7$Display, levels = rev(unique(TotalCounts$Display)))
  }
  if (!is.null(order_y_by)) {
    if ((length(order_y_by) == 1) && (order_y_by != "cluster")) {
      temp1 <- filter(abund7, .Group == order_y_by) %>%
        group_by(Display) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))

      abund7$Display <- factor(abund7$Display, levels = rev(temp1$Display))
    }
    if (length(order_y_by) > 1) {
      abund7$Display <- factor(abund7$Display, levels = order_y_by)
    }
    if ((length(order_y_by) == 1) && (order_y_by == "cluster")) {
      if (is.null(max_abundance)) {
        max_abundance <- max(abund7$Abundance)
      }
      tdata <- mutate(abund7,
        Abundance = ifelse(Abundance < min_abundance, min_abundance, Abundance),
        Abundance = ifelse(Abundance > max_abundance, max_abundance, Abundance)
      )
      tdata <- data.table::dcast(
        data.table::setDT(tdata),
        Display ~ .Group,
        value.var = "Abundance",
        fun.aggregate = sum
      )
      tdata <- as.matrix(tdata)
      rownames(tdata) <- tdata[, 1]
      tclust <- hclust(dist(tdata[, -1]))
      abund7$Display <- factor(abund7$Display, levels = tclust$labels[tclust$order])
    }
  }

  # show a row with the total sum of remaining taxa not among the top tax_show taxa
  if (isTRUE(showRemainingTaxa)) {
    uniqueTopTaxa <- unique(as.character(abund7$Display))
    nUniqueTopTaxa <- length(uniqueTopTaxa)
    nUniqueTotalTaxa <- length(unique(as.character(TotalCounts$Display)))

    if (nUniqueTopTaxa < nUniqueTotalTaxa) {
      remainingTaxa <- data.table::data.table(
        filter(
          abund6,
          !Display %in% unique(as.character(abund7$Display))
        )
      )
      remainingTaxa <- remainingTaxa[
        ,
        Display := paste0("Remaining taxa (", nUniqueTotalTaxa - nUniqueTopTaxa, ")")
      ][
        ,
        .(
          Abundance = sum(Abundance),
          Sum = sum(Sum),
          .Group = .Group[1]
        ),
        by = .(Display, Sample)
      ]
      abund7 <- as.data.frame(
        data.table::rbindlist(
          list(remainingTaxa, abund7)
        )
      )
    }
  }

  ## Order.x
  if (!is.null(order_x_by)) {
    if ((length(order_x_by) == 1) && (order_x_by != "cluster")) {
      temp1 <- filter(abund7, Display == order_x_by) %>%
        group_by(.Group) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      abund7$.Group <- factor(abund7$.Group, levels = as.character(temp1$.Group))
    }
    if (length(order_x_by) > 1) {
      abund7$.Group <- factor(abund7$.Group, levels = order_x_by)
    }
    if ((length(order_x_by) == 1) && (order_x_by == "cluster")) {
      if (is.null(max_abundance)) {
        max_abundance <- max(abund7$Abundance)
      }
      tdata <- mutate(abund7,
        Abundance = ifelse(Abundance < min_abundance, min_abundance, Abundance),
        Abundance = ifelse(Abundance > max_abundance, max_abundance, Abundance)
      )
      tdata <- data.table::dcast(data.table::setDT(tdata),
        Display ~ .Group,
        value.var = "Abundance",
        fun.aggregate = sum
      )
      tdata <- as.matrix(tdata)
      rownames(tdata) <- tdata[, 1]
      tclust <- hclust(dist(t(tdata[, -1])))
      abund7$.Group <- factor(abund7$.Group, levels = tclust$labels[tclust$order])
    }
  }

  ## Handle NA values
  if (plot_na == FALSE) {
    plot_na <- "grey50"
  } else {
    if (!is.null(color_vector)) {
      plot_na <- color_vector[1]
    } else {
      plot_na <- "#67A9CF"
    }
  }

  ## Scale to percentages if not normalised and scaled
  if (length(group_by) > 1) {
    abund7 <- merge(abund7, oldGroup, by = ".Group")
  }

  if (is.null(min_abundance)) {
    min_abundance <- ifelse(min(abund7$Abundance) > 0.001, min(abund7$Abundance), 0.001)
  }
  if (is.null(max_abundance)) {
    max_abundance <- max(abund7$Abundance)
  }

  ## Define the output
  if (!isTRUE(textmap)) {
    ## Make a heatmap style plot
    heatmap <- ggplot(abund7, aes_string(x = ".Group", y = "Display", label = formatC("Abundance", format = "f", digits = 1))) +
      geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) +
      theme(
        axis.text.y = element_text(size = 12, color = "black", vjust = 0.4),
        axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, angle = 90, hjust = 1),
        axis.title = element_blank(),
        text = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        # axis.ticks.length = unit(1, "mm"),
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        title = element_text(size = 8),
        panel.background = element_blank()
      )

    ## Get colorpalette for colorscale or set default
    if (!is.null(color_vector)) {
      color.pal <- color_vector
    } else {
      color.pal <- rev(RColorBrewer::brewer.pal(3, "RdYlBu"))
    }

    if (plot_values == TRUE) {
      abund8 <- abund7
      abund8$Abundance <- round(abund8$Abundance, round)
      heatmap <- heatmap + geom_text(data = abund8, size = plot_values_size, colour = "grey10", check_overlap = TRUE) +
        theme(legend.position = "none")
    }
    if (is.null(plot_legendbreaks)) {
      heatmap <- heatmap + scale_fill_gradientn(colours = color.pal, trans = plot_colorscale, na.value = plot_na, oob = scales::squish, limits = c(min_abundance, max_abundance))
    }
    if (!is.null(plot_legendbreaks)) {
      heatmap <- heatmap + scale_fill_gradientn(colours = color.pal, trans = plot_colorscale, breaks = plot_legendbreaks, na.value = plot_na, oob = scales::squish, limits = c(min_abundance, max_abundance))
    }


    if (is.null(normalise_by)) {
      heatmap <- heatmap + labs(x = "", y = "", fill = "% Read\nAbundance")
    }
    if (!is.null(normalise_by)) {
      heatmap <- heatmap + labs(x = "", y = "", fill = "Relative")
    }

    if (!is.null(facet_by)) {
      if (length(ogroup) > 1) {
        heatmap$data$.Group <- apply(heatmap$data[, ogroup], 1, paste, collapse = " ")
      } else {
        heatmap$data$.Group <- heatmap$data[, ogroup]
      }

      if (plot_values == TRUE) {
        if (length(ogroup) > 1) {
          heatmap$layers[[2]]$data$.Group <- apply(heatmap$layers[[2]]$data[, ogroup], 1, paste, collapse = " ")
        } else {
          heatmap$layers[[2]]$data$.Group <- heatmap$layers[[2]]$data[, ogroup]
        }
      }
      heatmap <- heatmap + facet_grid(reformulate(facet_by), scales = "free_x", space = "free")
      heatmap <- heatmap + theme(strip.text = element_text(size = 10))
    }

    # Return a function grid next to the heatmap with known functions about the Genera
    if (isTRUE(plot_functions)) {
      # Retrieve the genus names from the plot
      names <- data.frame(do.call("rbind", strsplit(levels(droplevels(heatmap$data$Display)), "; ", fixed = TRUE)))
      names <- data.frame(Genus = names[, which(c(tax_add, tax_aggregate) == "Genus")])
      names$Genus <- as.character(names$Genus)

      # Merge the genus and function information
      nameFunc <- merge(
        x = names,
        y = function_data[, c("Genus", functions)],
        all.x = TRUE,
        all.y = FALSE
      )
      nameFunc[is.na(nameFunc)] <- "NT"
      nameFuncM <- data.table::melt(
        data.table::as.data.table(nameFunc),
        id.vars = "Genus",
        value.name = "Value",
        variable.name = "Function"
      )
      nameFuncM$Value <- factor(nameFuncM$Value, levels = c("POS", "VAR", "NEG", "NT"))
      nameFuncM$Genus <- factor(nameFuncM$Genus, levels = names$Genus)

      functions_plot <- ggplot(nameFuncM, aes(x = Function, y = Genus, color = Value)) +
        geom_point(size = 4) +
        scale_color_manual(values = c("#31a354", "orange", "#f03b20", "grey90"), drop = FALSE) +
        theme(
          axis.text.x = element_text(
            size = 10,
            color = "black",
            angle = 90,
            hjust = 1,
            vjust = 0.4
          ),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.ticks.length = unit(1, "mm"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "grey95"),
          legend.key = element_blank()
        )
      # To be able to completely customise both heatmap and functions plot
      # individually by normal ggplot2 syntax (+ theme() etc...) before
      # creating the multiplot:
      amp_heatmap_function <- list(
        heatmap = heatmap,
        functions = functions_plot
      )
      structure(amp_heatmap_function,
        class = "hmfunplot",
        rel_widths = rel_widths
      )
    } else if (!isTRUE(plot_functions)) {
      return(heatmap)
    }
  } else if (isTRUE(textmap)) {
    # raw text heatmap data frame
    textmap <- abund7[, c("Display", "Abundance", ".Group"), drop = FALSE] %>%
      unique() %>%
      spread(key = .Group, value = Abundance)
    if (isTRUE(plot_functions)) {
      textmap <- merge(cbind(textmap,
        Genus = as.character(data.frame(do.call("rbind", strsplit(levels(droplevels(textmap$Display)), "; ", fixed = TRUE)))[, which(c(tax_add, tax_aggregate) == "Genus")])
      ),
      function_data[, c("Genus", switch(isTRUE(plot_functions),
        functions
      )), drop = FALSE],
      by = "Genus",
      all.x = TRUE,
      all.y = FALSE,
      fill = NA
      )
    }
    textmap <- textmap %>%
      arrange(desc(droplevels(Display)))
    textmap <- data.frame(textmap[, -c(which(colnames(textmap) %in% c("Display", "Genus"))), drop = FALSE],
      row.names = textmap$Display,
      check.names = FALSE
    )
    return(textmap)
  }
}
