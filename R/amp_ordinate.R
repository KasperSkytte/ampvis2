#' Ordination plot
#'
#' A wrapper around the vegan package to generate ggplot2 ordination plots suited for analysis and comparison of microbial communities. Simply choose an ordination type and a plot is returned.
#'
#' @usage amp_ordinate(data, type = "", transform = "", distmeasure = "", constrain = "")
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param filter_species Remove low abundant OTU's across all samples below this threshold in percent. Setting this to 0 may drastically increase computation time. (\emph{default}: \code{0.1})
#' @param type (\emph{required}) Type of ordination method. One of:
#' \itemize{
#'    \item \code{"PCA"}: (\emph{default}) Principal Components Analysis
#'    \item \code{"RDA"}: Redundancy Analysis (considered the constrained version of PCA)
#'    \item \code{"CA"}: Correspondence Analysis
#'    \item \code{"CCA"}: Canonical Correspondence Analysis (considered the constrained version of CA)
#'    \item \code{"DCA"}: Detrended Correspondence Analysis
#'    \item \code{"NMDS"}: non-metric Multidimensional Scaling
#'    \item \code{"PCOA"} or \code{"MMDS"}: metric Multidimensional Scaling a.k.a Principal Coordinates Analysis (not to be confused with PCA)
#'    }
#'    \emph{Note that PCoA is not performed by the vegan package, but the \code{\link[ape]{pcoa}} function from the APE package.}
#' @param distmeasure (\emph{required for nMDS and PCoA}) Distance measure used for the distance-based ordination methods (nMDS and PCoA). Choose one of the following:
#' \itemize{
#'   \item \code{"wunifrac"}: Weighted generalized UniFrac distances (alpha=0.5), calculated by \code{\link[GUniFrac]{GUniFrac}}. Requires a phylogenetic tree.
#'   \item \code{"unifrac"}: Unweighted UniFrac distances, calculated by \code{\link[GUniFrac]{GUniFrac}}. Requires a phylogenetic tree.
#'   \item \code{"jsd"}: Jensen-Shannon Divergence, based on \url{http://enterotype.embl.de/enterotypes.html}.
#'   \item Any of the distance measures supported by \code{\link[vegan]{vegdist}}: \code{"manhattan"}, \code{"euclidean"}, \code{"canberra"}, \code{"bray"}, \code{"kulczynski"}, \code{"jaccard"}, \code{"gower"}, \code{"altGower"}, \code{"morisita"}, \code{"horn"}, \code{"mountford"}, \code{"raup"}, \code{"binomial"}, \code{"chao"}, \code{"cao"}, \code{"mahalanobis"}.
#'   \item or \code{"none"}. (\emph{default})
#'  }
#' You can also write your own math formula, see details in \code{\link[vegan]{vegdist}}.
#' @param transform (\emph{recommended}) Transforms the abundances before ordination, choose one of the following: \code{"total"}, \code{"max"}, \code{"freq"}, \code{"normalize"}, \code{"range"}, \code{"standardize"}, \code{"pa"} (presence/absense), \code{"chi.square"}, \code{"hellinger"}, \code{"log"}, or \code{"sqrt"}, see details in \code{\link[vegan]{decostand}}. Using the hellinger transformation is a good choice when performing PCA/RDA as it will produce a more ecologically meaningful result (read about the double-zero problem in Numerical Ecology). When the Hellinger transformation is used with CA/CCA it will help reducing the impact of low abundant species. When performing nMDS or PCoA (aka mMDS) it is not recommended to also use data transformation as this will obscure the chosen distance measure. (\emph{default:} \code{"hellinger"})
#' @param constrain (\emph{required for RDA and CCA}) Variable(s) in the metadata for constrained analyses (RDA and CCA). Multiple variables can be provided by a vector, fx \code{c("Year", "Temperature")}, but keep in mind that the more variables selected the more the result will be similar to unconstrained analysis.
#' @param x_axis Which axis from the ordination results to plot as the first axis. Have a look at the \code{$screeplot} with \code{detailed_output = TRUE} to validate axes. (\emph{default:} \code{1})
#' @param y_axis Which axis from the ordination results to plot as the second axis. Have a look at the \code{$screeplot} with \code{detailed_output = TRUE} to validate axes. (\emph{default:} \code{2})
#' @param print_caption Auto-generate a figure caption based on the arguments used. The caption includes a description of how the result has been generated as well as references for the methods used.
#' @param sample_color_by Color sample points by a variable in the metadata.
#' @param sample_color_order Order the colors in \code{sample_color_by} by the order in a vector.
#' @param sample_label_by Label sample points by a variable in the metadata.
#' @param sample_label_size Sample labels text size. (\emph{default:} \code{4})
#' @param sample_label_segment_color Sample labels repel-segment color. (\emph{default:} \code{"black"})
#' @param sample_shape_by Shape sample points by a variable in the metadata.
#' @param sample_colorframe Frame the sample points with a polygon by a variable in the metadata split by the variable defined by \code{sample_color_by}, or simply \code{TRUE} to frame the points colored by \code{sample_color_by}. (\emph{default:} \code{FALSE})
#' @param sample_colorframe_label Label by a variable in the metadata.
#' @param sample_point_size Size of the sample points. (\emph{default:} \code{2})
#' @param sample_trajectory Make a trajectory between sample points by a variable in the metadata.
#' @param sample_trajectory_group Make a trajectory between sample points by the \code{sample_trajectory} argument, but within individual groups.
#' @param sample_plotly Enable interactive sample points so that they can be hovered to show additional information from the metadata. Provide a vector of the metadata variables to show, or \code{"all"} to display all. Click or double click the elements in the legend to hide/show parts of the data. To hide the legend use \code{plotly::layout(amp_ordinate(...), showlegend = FALSE)}, see more options at \url{https://plot.ly/r/}.
#'
#' @param species_plot (\emph{logical}) Plot species points or not. (\emph{default:} \code{FALSE})
#' @param species_shape The shape of the species points, fx \code{1} for hollow circles or \code{20} for dots. (\emph{default:} \code{20})
#' @param species_point_size Size of the species points. (\emph{default:} \code{2})
#' @param species_nlabels Number of the most extreme species labels to plot (ordered by the sum of the numerical values of the x,y coordinates. Only makes sense with PCA/RDA).
#' @param species_label_taxonomy Taxonomic level by which to label the species points. (\emph{default:} \code{"Genus"})
#' @param species_label_size Size of the species text labels. (\emph{default:} \code{3})
#' @param species_label_color Color of the species text labels. (\emph{default:} \code{"grey10"})
#' @param species_rescale (\emph{logical}) Rescale species points or not. Basically they will be multiplied by 0.8, for visual convenience only. (\emph{default:} \code{FALSE})
#' @param species_plotly (\emph{logical}) Enable interactive species points so that they can be hovered to show complete taxonomic information. (\emph{default:} \code{FALSE})
#'
#' @param envfit_factor A vector of categorical environmental variables from the metadata to fit onto the ordination plot. See details in \code{\link[vegan]{envfit}}.
#' @param envfit_numeric A vector of numerical environmental variables from the metadata to fit arrows onto the ordination plot. The lengths of the arrows are scaled by significance. See details in \code{\link[vegan]{envfit}}.
#' @param envfit_signif_level The significance threshold for displaying the results of \code{envfit_factor} or \code{envfit_numeric}. (\emph{default:} \code{0.005})
#' @param envfit_textsize Size of the envfit text on the plot. (\emph{default:} \code{3})
#' @param envfit_textcolor Color of the envfit text on the plot. (\emph{default:} \code{"darkred"})
#' @param envfit_numeric_arrows_scale Scale the size of the numeric arrows. (\emph{default:} \code{1})
#' @param envfit_arrowcolor Color of the envfit arrows on the plot. (\emph{default:} \code{"darkred"})
#' @param envfit_show (\emph{logical}) Show the results on the plot or not. (\emph{default:} \code{TRUE})
#'
#' @param repel_labels (\emph{logical}) Repel all labels to prevent cluttering of the plot. (\emph{default:} \code{TRUE})
#' @param opacity Opacity of all plotted points and sample_colorframe. \code{0}: invisible, \code{1}: opaque. (\emph{default:} \code{0.8})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param detailed_output (\emph{logical}) Return additional details or not (model, scores, figure caption, inputmatrix, screeplot etc). If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' @param ... Pass additional arguments to the vegan ordination functions, fx the \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{metaMDS}} functions, see the documentation.
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @details
#' The \code{\code{amp_ordinate}} function is primarily based on two packages; \code{\link[vegan]{vegan-package}}, which performs the actual ordination, and the \code{\link[ggplot2]{ggplot2-package}} to generate the plot. The function generates an ordination plot by the following process:
#' \enumerate{
#'   \item Various input argument checks and error messages
#'   \item OTU-table filtering, where low abundant OTU's across all samples are removed (if not \code{filter_species = 0} is set)
#'   \item Data transformation (if not \code{transform = "none"} is set)
#'   \item Calculate distance matrix based on the chosen \code{distmeasure} if the chosen ordination method is PCoA/nMDS/DCA
#'   \item Perform the actual ordination and calculate the axis scores for both samples and species/OTU's
#'   \item Visualise the result with ggplot2 or plotly in various ways defined by the user
#' }
#'
#'
#' @section Using a custom distance matrix:
#' If you wan't to calculate a distance matrix manually and use it for PCoA or nMDS in \code{amp_ordinate}, it can be done by setting \code{filter_species = 0}, \code{transform = "none"}, \code{distmeasure = "none"}, and then override the abundance table (\code{$abund}) in the ampvis2 object, like below. The matrix must be a symmetrical matrix containing coefficients for all pairs of samples in the data.
#' \preformatted{
#' #Override the abundance table in the ampvis2 object with a custom distance matrix
#' ampvis2_object$abund <- custom_dist_matrix
#' #set filter_species = 0, transform = "none", and distmeasure = "none"
#' amp_ordinate(ampvis2_object,
#'              type = "pcoa",
#'              filter_species = 0,
#'              transform = "none",
#'              distmeasure = "none")
#' }
#'
#' @import ggplot2
#' @importFrom dplyr arrange group_by summarise desc
#' @importFrom plotly ggplotly
#' @importFrom GUniFrac GUniFrac
#' @importFrom ggrepel geom_text_repel
#' @importFrom ape pcoa
#' @importFrom vegan cca decorana decostand envfit metaMDS rda scores vegdist
#' @importFrom purrr imap
#' @importFrom tidyr unite
#'
#' @export
#'
#' @references
#'   GUide to STatistical Analysis in Microbial Ecology (GUSTA ME): \url{https://mb3is.megx.net/gustame}
#'
#'   Legendre, Pierre & Legendre, Louis (2012). Numerical Ecology. Elsevier Science. ISBN: 9780444538680
#'
#'   Legendre, P., & Gallagher, E. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia, 129(2), 271-280. \url{http://doi.org/10.1007/s004420100716}
#'
#' @seealso
#' \code{\link{amp_load}}
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#' 
#' # PCA with data transformation, colored by WWTP
#' amp_ordinate(AalborgWWTPs,
#'   type = "PCA",
#'   transform = "hellinger",
#'   sample_color_by = "Plant",
#'   sample_colorframe = TRUE
#' )
#' \dontrun{
#' # Interactive CCA with data transformation constrained to seasonal period
#' amp_ordinate(AalborgWWTPs,
#'   type = "CCA",
#'   transform = "Hellinger",
#'   constrain = "Period",
#'   sample_color_by = "Period",
#'   sample_colorframe = TRUE,
#'   sample_colorframe_label = "Period",
#'   sample_plotly = "all"
#' )
#' }
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_ordinate <- function(data,
                         filter_species = 0.1,
                         type = "PCA",
                         distmeasure = "none",
                         transform = "hellinger",
                         constrain = NULL,
                         x_axis = 1,
                         y_axis = 2,
                         print_caption = FALSE,
                         sample_color_by = NULL,
                         sample_color_order = NULL,
                         sample_shape_by = NULL,
                         sample_colorframe = FALSE,
                         sample_colorframe_label = NULL,
                         sample_label_by = NULL,
                         sample_label_size = 4,
                         sample_label_segment_color = "black",
                         sample_point_size = 2,
                         sample_trajectory = NULL,
                         sample_trajectory_group = NULL,
                         sample_plotly = NULL,
                         species_plot = FALSE,
                         species_nlabels = 0,
                         species_label_taxonomy = "Genus",
                         species_label_size = 3,
                         species_label_color = "grey10",
                         species_rescale = FALSE,
                         species_point_size = 2,
                         species_shape = 20,
                         species_plotly = FALSE,
                         envfit_factor = NULL,
                         envfit_numeric = NULL,
                         envfit_signif_level = 0.005,
                         envfit_textsize = 3,
                         envfit_textcolor = "darkred",
                         envfit_numeric_arrows_scale = 1,
                         envfit_arrowcolor = "darkred",
                         envfit_show = TRUE,
                         repel_labels = TRUE,
                         opacity = 0.8,
                         tax_empty = "best",
                         detailed_output = FALSE,
                         ...) {

  ### Data must be in ampvis2 format
  if (class(data) != "ampvis2") {
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  }

  ##### Sanity check of options  #####
  if (species_plotly == TRUE & !is.null(sample_plotly)) {
    stop("You can not use plotly for both species and samples in the same plot.", call. = FALSE)
  }
  if (species_plotly == TRUE | !is.null(sample_plotly)) {
    # message("geom_text_repel is not supported by plotly yet.")
    repel_labels <- FALSE
  }

  # Impossible to do ordination with 1 or 2 samples
  if (length(unique(data$metadata[, 1])) <= 2) {
    stop("Ordination cannot be performed on 2 or fewer samples (the number of resulting axes will always be n-1, where n is the number of samples).", call. = FALSE)
  }

  if (is.null(sample_color_by) & !is.logical(sample_colorframe) & !is.null(sample_colorframe)) {
    sample_color_by <- sample_colorframe
  }

  # Check the data
  data <- amp_rename(data = data, tax_empty = tax_empty)

  ##### Filter #####
  if (filter_species != 0 && is.numeric(filter_species)) {
    # First transform to percentages
    abund_pct <- data$abund
    abund_pct[, which(colSums(abund_pct) != 0)] <- as.data.frame(apply(abund_pct[, which(colSums(abund_pct) != 0), drop = FALSE], 2, function(x) x / sum(x) * 100))
    rownames(abund_pct) <- rownames(data$abund) # keep rownames

    # Then filter low abundant OTU's where ALL samples have below the threshold set with filter_species in percent
    abund_subset <- abund_pct[!apply(abund_pct, 1, function(row) all(row <= filter_species)), , drop = FALSE] # remove low abundant OTU's
    data$abund <- data$abund[which(rownames(data$abund) %in% rownames(abund_subset)), , drop = FALSE]
    rownames(data$tax) <- data$tax$OTU
    data$tax <- data$tax[which(rownames(data$tax) %in% rownames(abund_subset)), , drop = FALSE] # same with taxonomy
  }

  # to fix user argument characters, so fx PCoA/PCOA/pcoa are all valid
  type <- tolower(type)
  transform <- tolower(transform)
  distmeasure <- tolower(distmeasure)

  ##### Data transformation with decostand()  #####
  if (!transform == "none" & transform != "sqrt") {
    transform <- tolower(transform)
    data$abund <- t(vegan::decostand(t(data$abund), method = transform))
  } else if (transform == "sqrt") {
    data$abund <- t(sqrt(t(data$abund)))
  }

  ##### Inputmatrix AFTER transformation  #####
  if (any(type == c("nmds", "mmds", "pcoa", "dca"))) {
    if (!type == "nmds" & (species_plot == TRUE | species_plotly == TRUE)) {
      stop("No speciesscores available with mMDS/PCoA, DCA.", call. = FALSE)
    }
    if (!distmeasure == "none") {
      # Calculate distance matrix with vegdist()
      if (distmeasure == "jsd") {
        # This is based on http://enterotype.embl.de/enterotypes.html
        # Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
        # Unfortunately this code is SLOOOOOOOOW
        dist.JSD <- function(inMatrix, pseudocount = 0.000001) {
          KLD <- function(x, y) sum(x * log(x / y))
          JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
          matrixColSize <- length(colnames(inMatrix))
          matrixRowSize <- length(rownames(inMatrix))
          colnames <- colnames(inMatrix)
          resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

          inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))

          for (i in 1:matrixColSize) {
            for (j in 1:matrixColSize) {
              resultsMatrix[i, j] <- JSD(
                as.vector(inMatrix[, i]),
                as.vector(inMatrix[, j])
              )
            }
          }
          colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
          as.dist(resultsMatrix) -> resultsMatrix
          attr(resultsMatrix, "method") <- "dist"
          return(resultsMatrix)
        }
        message("Calculating Jensen-Shannon Divergence (JSD) distances... ")
        inputmatrix <- dist.JSD(data$abund)
        message("Done.")
      } else if (any(distmeasure == c(
        "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower",
        "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"
      ))) {
        message("Calculating distance matrix... ")
        inputmatrix <- vegan::vegdist(t(data$abund), method = distmeasure)
        message("Done.")
      } else if (distmeasure == "unifrac") {
        # no tree no unifrac
        if (!any(names(data) == "tree")) {
          stop("No phylogenetic tree in the provided data.", call. = FALSE)
        }
        message("Calculating unweighted generalized UniFrac distances... ")
        unifracs <- suppressWarnings(GUniFrac::GUniFrac(t(data$abund), data$tree, alpha = 0.5)$unifracs)
        message("Done.")
        # GUniFrac returns both weighted and unweighted distances in the same dataframe
        inputmatrix <- unifracs[, , "d_UW"] %>%
          as.data.frame()
      } else if (distmeasure == "wunifrac") {
        # no tree no unifrac
        if (!any(names(data) == "tree")) {
          stop("No phylogenetic tree in the provided data.", call. = FALSE)
        }
        message("Calculating weighted generalized UniFrac distances (alpha=0.5)... ")
        unifracs <- suppressWarnings(GUniFrac::GUniFrac(t(data$abund), data$tree, alpha = 0.5)$unifracs)
        message("Done.")
        # GUniFrac returns both weighted and unweighted distances in the same dataframe
        inputmatrix <- unifracs[, , "d_0.5"] %>%
          as.data.frame()
      }
    } else if (distmeasure == "none") {
      warning("No distance measure selected, using raw data. If this is not deliberate, please provide one with the argument: distmeasure.", call. = FALSE)
      inputmatrix <- t(data$abund)
    }

    if (transform != "none" & distmeasure != "none") {
      warning("Using both transformation AND a distance measure is not recommended for distance-based ordination (nMDS/PCoA/DCA). If this is not deliberate, consider transform = \"none\".", call. = FALSE)
    }
  } else if (any(type == c("pca", "rda", "ca", "cca"))) {
    inputmatrix <- t(data$abund)
  }

  ##### Perform ordination  #####
  # Generate data depending on the chosen ordination type
  if (type == "pca") {
    # make the model
    model <- vegan::rda(inputmatrix, ...)

    # axis (and data column) names
    x_axis_name <- paste0("PC", x_axis)
    y_axis_name <- paste0("PC", y_axis)

    # Calculate the amount of inertia explained by each axis
    totalvar <- round(model$CA$eig / model$tot.chi * 100, 1)

    # Calculate species- and site scores
    sitescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if (type == "rda") {
    if (is.null(constrain)) {
      stop("Argument constrain must be provided when performing constrained/canonical analysis.", call. = FALSE)
    }
    # make the model
    codestring <- paste0("rda(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") # function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <- eval(parse(text = codestring))

    # axes depend on the results
    x_axis_name <- paste0("RDA", x_axis)
    if (model$CCA$rank <= 1) {
      y_axis_name <- "PC1"

      # Calculate the amount of inertia explained by each axis
      totalvar <- c(
        round(model$CCA$eig / model$tot.chi * 100, 1), # constrained of total space
        round(model$CA$eig / model$tot.chi * 100, 1) # UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig / model$CA$tot.chi * 100, 1)) # UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("RDA", y_axis)

      # Calculate the amount of inertia explained by each axis
      totalvar <- c(
        round(model$CCA$eig / model$tot.chi * 100, 1), # constrained of total space
        round(model$CA$eig / model$tot.chi * 100, 1) # UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig / model$CCA$tot.chi * 100, 1)) # constrained of total constrained space
    }

    # Calculate species- and site scores
    sitescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if (type == "nmds") {
    # make the model
    if (ncol(data$abund) > 100) {
      message("Performing non-Metric Multidimensional Scaling on more than 100 samples, this may take some time ... ")
    }
    model <- vegan::metaMDS(inputmatrix, trace = FALSE, ...)
    if (ncol(data$abund) > 100) {
      message("Done.")
    }

    # axis (and data column) names
    x_axis_name <- paste0("NMDS", x_axis)
    y_axis_name <- paste0("NMDS", y_axis)

    # Calculate species- and site scores
    # Speciesscores may not be available with MDS
    sitescores <- vegan::scores(model, display = "sites")
    if (!length(model$species) > 1) {
      speciesscores <- NULL
      if (species_plot == TRUE | species_plotly == TRUE) {
        stop("Speciesscores are not available.", call. = FALSE)
      }
    } else {
      speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
    }
  } else if (type == "mmds" | type == "pcoa") {
    # make the model
    model <- ape::pcoa(inputmatrix, ...)

    # axis (and data column) names
    x_axis_name <- paste0("PCo", x_axis)
    y_axis_name <- paste0("PCo", y_axis)

    # Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$values$Relative_eig * 100, 1)
    names(totalvar) <- c(paste0("PCo", seq(1:length(totalvar))))

    # Calculate species- and site scores
    # Speciesscores are not available with pcoa
    sitescores <- as.data.frame(model$vectors)
    colnames(sitescores) <- c(paste0("PCo", seq(1:length(sitescores))))
    speciesscores <- NULL
  } else if (type == "ca") {
    # make the model
    model <- vegan::cca(inputmatrix, ...)

    # axis (and data column) names
    x_axis_name <- paste0("CA", x_axis)
    y_axis_name <- paste0("CA", y_axis)

    # Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$CA$eig / model$tot.chi * 100, 1)

    # Calculate species- and site scores
    sitescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if (type == "cca") {
    if (is.null(constrain)) {
      stop("Argument constrain must be provided when performing constrained/canonical analysis.", call. = FALSE)
    }
    # make the model
    codestring <- paste0("cca(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") # function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <- eval(parse(text = codestring))

    # axes depend on the results
    x_axis_name <- paste0("CCA", x_axis)
    if (model$CCA$rank <= 1) {
      y_axis_name <- "CA1"

      # Calculate the amount of inertia explained by each axis
      totalvar <- c(
        round(model$CCA$eig / model$tot.chi * 100, 1), # constrained of total space
        round(model$CA$eig / model$tot.chi * 100, 1) # UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig / model$CA$tot.chi * 100, 1)) # UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("CCA", y_axis)

      # Calculate the amount of inertia explained by each axis
      totalvar <- c(
        round(model$CCA$eig / model$tot.chi * 100, 1), # constrained of total space
        round(model$CA$eig / model$tot.chi * 100, 1) # UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig / model$CCA$tot.chi * 100, 1)) # constrained of total constrained space
    }

    # Calculate species- and site scores
    sitescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if (type == "dca") {
    # make the model
    model <- vegan::decorana(inputmatrix, ...)

    # axis (and data column) names
    x_axis_name <- paste0("DCA", x_axis)
    y_axis_name <- paste0("DCA", y_axis)

    # Calculate the percentage of eigenvalues explained by the axes
    # totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)

    # Calculate species- and site scores
    sitescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  }

  ##### Data for ggplot  #####
  dsites <- cbind.data.frame(data$metadata, sitescores)

  if (!is.null(sample_color_order)) {
    dsites[, sample_color_by] <- factor(dsites[, sample_color_by], levels = sample_color_order)
  }

  if (length(speciesscores) > 1) {
    dspecies <- merge(data.frame(speciesscores, OTU = rownames(speciesscores)), data$tax, by.x = "OTU")
    dspecies$dist <- dspecies[, x_axis_name]^2 + dspecies[, y_axis_name]^2
    dspecies <- dplyr::arrange(dspecies, desc(dist))
    rownames(dspecies) <- dspecies$OTU
    if (species_rescale == TRUE) {
      maxx <- max(abs(dsites[, x_axis_name])) / max(abs(dspecies[, x_axis_name]))
      dspecies[, x_axis_name] <- dspecies[, x_axis_name] * maxx * 0.8
      maxy <- max(abs(dsites[, y_axis_name])) / max(abs(dspecies[, y_axis_name]))
      dspecies[, y_axis_name] <- dspecies[, y_axis_name] * maxy * 0.8
    }
  } else {
    dspecies <- NULL
  }

  ##### Base plot object #####
  if (is.null(sample_color_by) && is.null(sample_shape_by)) {
    plot <- ggplot(
      dsites,
      aes_string(
        x = x_axis_name,
        y = y_axis_name
      )
    )
  } else if (!is.null(sample_color_by) && is.null(sample_shape_by)) {
    plot <- ggplot(
      dsites,
      aes_string(
        x = x_axis_name,
        y = y_axis_name,
        color = sample_color_by
      )
    )
  } else if (is.null(sample_color_by) && !is.null(sample_shape_by)) {
    plot <- ggplot(
      dsites,
      aes_string(
        x = x_axis_name,
        y = y_axis_name,
        shape = sample_shape_by
      )
    )
  } else if (!is.null(sample_color_by) && !is.null(sample_shape_by)) {
    plot <- ggplot(
      dsites,
      aes_string(
        x = x_axis_name,
        y = y_axis_name,
        color = sample_color_by,
        shape = sample_shape_by
      )
    )
  }

  ##### Colorframe  #####
  if (!sample_colorframe == FALSE) {
    if (is.null(sample_color_by) & sample_colorframe == TRUE) {
      stop("Applying a colorframe to the sample points requires a variable in the metadata to be provided by the argument sample_colorframe, or by sample_color_by if the former is only TRUE.", call. = FALSE)
    }
    if (sample_colorframe == TRUE) {
      splitData <- base::split(plot$data, plot$data[, sample_color_by]) %>%
        lapply(function(df) {
          df[chull(df[, x_axis_name], df[, y_axis_name]), ]
        })
      hulls <- do.call(rbind, splitData)
      plot <- plot + geom_polygon(data = hulls, aes_string(fill = sample_color_by, group = sample_color_by), alpha = 0.2 * opacity)
    } else if (!is.logical(sample_colorframe) & !is.null(sample_colorframe)) {
      plot$data$colorframeGroup <- paste(plot$data[, sample_color_by], plot$data[, sample_colorframe]) %>% as.factor()
      splitData <- base::split(plot$data, plot$data$colorframeGroup) %>%
        lapply(function(df) {
          df[chull(df[, x_axis_name], df[, y_axis_name]), ]
        })
      hulls <- do.call(rbind, splitData)
      plot <- plot + geom_polygon(data = hulls, aes_string(fill = sample_color_by, group = "colorframeGroup"), alpha = 0.2 * opacity)
    }
  }

  ##### Plot sample points  #####
  if (!is.null(sample_plotly)) {
    if (identical(sample_plotly, "all") | isTRUE(sample_plotly)) {
      sample_plotly <- colnames(data$metadata)
    }
    data_plotly <- apply(
      data$metadata[, sample_plotly, drop = FALSE],
      1,
      function(x) {
        paste0(paste0(names(x), ": "), x, collapse = "<br>")
      }
    )
    plot <- plot +
      suppressWarnings(geom_point(
        size = 2, alpha = opacity,
        aes(text = data_plotly)
      )) + # HER
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  } else {
    plot <- plot +
      geom_point(size = sample_point_size, alpha = opacity) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  }

  # Only eigenvalue-based ordination methods can be displayed with % on axes
  if (type == "pca" | type == "ca" | type == "pcoa" | type == "mmds") {
    plot <- plot +
      xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
      ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
  } else if (type == "rda" | type == "cca") {
    if (model$CCA$rank > 1) {
      plot <- plot +
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "% / ", constrainedvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    } else if (model$CCA$rank <= 1) {
      plot <- plot +
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    }
  } else if (type == "nmds") {
    plot <- plot +
      annotate("text", size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste0("Stress value = ", round(model$stress, 3)))
  }

  ##### Plot species points  #####
  if (species_plot == TRUE) {
    if (species_plotly == T) {
      # generate hover text for OTU's
      data_plotly <- data$tax %>%
        purrr::imap(~paste(.y, .x, sep = ": ")) %>%
        as.data.frame() %>%
        tidyr::unite("test", sep = "<br>") %>%
        unlist(use.names = FALSE) %>%
        unique()
      plot <- plot +
        geom_point(
          data = dspecies,
          color = "darkgrey",
          shape = species_shape,
          size = species_point_size - 1,
          alpha = opacity,
          aes(text = data_plotly)
        )
    } else {
      plot <- plot +
        geom_point(
          data = dspecies,
          color = "darkgrey",
          shape = species_shape,
          size = species_point_size,
          alpha = opacity
        )
    }
  }

  ##### Plot text labels  #####
  if (!is.null(sample_colorframe_label)) {
    temp <- data.frame(
      group = dsites[, sample_colorframe_label],
      x = dsites[, x_axis_name],
      y = dsites[, y_axis_name]
    ) %>%
      group_by(group) %>%
      summarise(cx = mean(x), cy = mean(y)) %>%
      as.data.frame()
    temp2 <- merge(dsites, temp,
      by.x = sample_colorframe_label,
      by.y = "group"
    )
    temp3 <- temp2[!duplicated(temp2[, sample_colorframe_label]), ]
    if (repel_labels == T) {
      plot <- plot + ggrepel::geom_text_repel(
        data = temp3,
        aes_string(
          x = "cx",
          y = "cy",
          label = sample_colorframe_label
        ),
        size = 3,
        color = "black",
        fontface = 2
      )
    } else {
      plot <- plot + geom_text(
        data = temp3,
        aes_string(
          x = "cx",
          y = "cy",
          label = sample_colorframe_label
        ),
        size = 3,
        color = "black",
        fontface = 2
      )
    }
  }

  ##### Sample_trajectory  #####
  if (!is.null(sample_trajectory)) {
    traj <- dsites[order(dsites[, sample_trajectory, drop = FALSE]), ]
    if (is.null(sample_trajectory_group)) {
      plot <- plot + geom_path(data = traj)
    } else if (!is.null(sample_trajectory_group)) {
      plot <- plot + geom_path(data = traj, aes_string(group = sample_trajectory_group))
    }
  }

  ##### Sample point labels  #####
  if (!is.null(sample_label_by)) {
    if (repel_labels == T) {
      plot <- plot + ggrepel::geom_text_repel(aes_string(label = sample_label_by), size = sample_label_size, color = "grey40", segment.color = sample_label_segment_color)
    }
    else {
      plot <- plot + geom_text(aes_string(label = sample_label_by), size = sample_label_size, color = "grey40", segment.color = sample_label_segment_color)
    }
  }

  ##### Plot species labels  #####
  if (species_nlabels > 0) {
    if (repel_labels == T) {
      plot <- plot + ggrepel::geom_text_repel(data = dspecies[1:species_nlabels, ], aes_string(x = x_axis_name, y = y_axis_name, label = species_label_taxonomy), colour = species_label_color, size = species_label_size, fontface = 4, inherit.aes = FALSE)
    }
    else {
      plot <- plot + geom_text(data = dspecies[1:species_nlabels, ], aes_string(x = x_axis_name, y = y_axis_name, label = species_label_taxonomy), colour = species_label_color, size = species_label_size, fontface = 4, inherit.aes = FALSE)
    }
  }

  ##### Categorical fitting  #####
  if (!is.null(envfit_factor)) {
    evf_factor_model <- envfit(dsites[, c(x_axis_name, y_axis_name)],
      data$metadata[, envfit_factor, drop = FALSE],
      permutations = 999
    )
    evf_factor_data <- data.frame(
      Name = rownames(evf_factor_model$factors$centroids),
      Variable = evf_factor_model$factors$var.id,
      evf_factor_model$factors$centroids,
      pval = evf_factor_model$factors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_factor_data) > 0 & envfit_show == TRUE) {
      if (repel_labels == T) {
        plot <- plot + ggrepel::geom_text_repel(data = evf_factor_data, aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_textcolor, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")
      }
      else {
        plot <- plot + geom_text(data = evf_factor_data, aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_textcolor, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")
      }
    }
    if (nrow(evf_factor_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.", call. = FALSE)
    }
  } else {
    evf_factor_model <- NULL
  }

  ##### Numerical fitting  #####
  if (!is.null(envfit_numeric)) {
    evf_numeric_model <- envfit(dsites[, c(x_axis_name, y_axis_name)],
      data$metadata[, envfit_numeric, drop = FALSE],
      permutations = 999
    )
    evf_numeric_data <- data.frame(
      Name = rownames(evf_numeric_model$vectors$arrows),
      evf_numeric_model$vectors$arrows * sqrt(evf_numeric_model$vectors$r) * envfit_numeric_arrows_scale,
      pval = evf_numeric_model$vectors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_numeric_data) > 0 & envfit_show == TRUE) {
      plot <- plot + geom_segment(
        data = evf_numeric_data,
        aes_string(
          x = 0,
          xend = x_axis_name,
          y = 0,
          yend = y_axis_name
        ),
        arrow = arrow(length = unit(3, "mm")),
        colour = envfit_arrowcolor,
        size = 1,
        inherit.aes = FALSE
      ) +
        geom_text(
          data = evf_numeric_data,
          aes_string(
            x = x_axis_name,
            y = y_axis_name,
            label = "Name"
          ),
          colour = envfit_textcolor,
          inherit.aes = FALSE,
          size = envfit_textsize,
          hjust = 1.2,
          vjust = 1.2,
          fontface = "bold"
        )
    }
    if (nrow(evf_numeric_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.", call. = FALSE)
    }
  } else {
    evf_numeric_model <- NULL
  }

  ##### generate figure caption #####
  # Generate caption
  caption <- paste0(
    switch(type,
      "pca" = "Principal Components Analysis (PCA)",
      "rda" = "Redundancy Analysis (RDA)",
      "ca" = "Correspondence Analysis (CA)",
      "cca" = "Canonical Correspondence Analysis (CCA)",
      "dca" = "Detrended Correspondence Analysis (DCA)",
      "nmds" = "non-Metric Multidimensional Scaling Analysis (nMDS)",
      "pcoa" = "Principal Coordinates Analysis (PCoA)",
      "mmds" = "metric Multidimensional Scaling (aka Principal Coordinates Analysis, PCoA)"
    ),
    if (any(type %in% c("nmds", "mmds", "pcoa", "dca"))) {
      paste0(
        " based on the ",
        switch(distmeasure,
          "manhattan" = "Manhattan (city-block) distance measure",
          "euclidean" = "Euclidean distance measure",
          "canberra" = "Canberra distance measure (Lance & Williams, 1966)",
          "clark" = "Clark distance measure (Clark & Evans, 1954)",
          "bray" = "Bray-Curtis distance measure (Bray & Curtis, 1957)",
          "kulczynski" = "Kulczynski distance measure (Kulczynski, 1928)",
          "jaccard" = "Jaccard distance measure (Jaccard, 1901)",
          "gower" = "Gower distance measure (Gower, 1971)",
          "altGower" = "modified Gower distance measure (Anderson et al., 2006)",
          "morisita" = "morisita distance measure (Morisita, 1959)",
          "horn" = "Horn-Morisita distance measure (Horn, 1966)",
          "mountford" = "Mountford distance measure (square-root of that described in Mountford, 1962)",
          "raup" = "Raup-Crick distance measure (Raup & Crick, 1979)",
          "binomial" = "Binomial distance measure (Anderson & Millar, 2004)",
          "chao" = "Chao distance measure (Chao et al., 2005)",
          "cao" = "Cao (CYd) distance measure (Cao et al, 1997)",
          "mahalanobis" = "Mahalanobis distance measure (Mahalanobis, 1936)",
          "none" = " using no distance measure"
        )
      )
    },
    " of ", ncol(data$abund), " samples",
    if (any(type %in% c("cca", "rda"))) {
      paste0(
        " constrained to the ",
        if (length(constrain) == 1) {
          paste0("\"", constrain, "\" variable")
        } else if (length(constrain) > 1) {
          paste0("variables: ", paste0(constrain, collapse = "+"))
        }
      )
    },
    ".",
    if (filter_species != 0 && is.numeric(filter_species)) {
      paste0(
        " Prior to the analysis, OTU's that are not present in more than ",
        filter_species, "% relative abundance in any sample have been removed."
      )
    },
    if (transform != "none" && any(transform %in% c(
      "total", "max", "freq", "normalize", "range",
      "rank", "rrank", "standardize", "pa", "chi.square",
      "hellinger", "log", "sqrt"
    ))) {
      switch(transform,
        "total" = " The data has been transformed initially by dividing the OTU counts in each sample by the total counts of all OTU's in the particular sample.",
        "max" = " The data has been transformed initially by dividing each OTU count by the maximum of the particular OTU in any sample.",
        "freq" = " The data has been transformed initially by dividing each OTU count by the maximum of the particular OTU in any sample, and then multiplied by the number of non-zero OTU counts, so that the average of non-zero OTU's is 1 (Oksanen, 1983).",
        "normalize" = " The data has been normalized initially by making the sum-of-squares of the OTU counts in each sample equal to 1.",
        "range" = " The data has been normalized initially to range from 0 to 1.",
        "rank" = " The data has been transformed initially by replacing the OTU counts with their increasing ranks, where zeros are unchanged and the average ranks are used for tied values.",
        "rrank" = " The data has been transformed initially by replacing the OTU counts with their increasing relative ranks (with max value 1), where zeros are unchanged and the average ranks are used for tied values.",
        "standardize" = " The data has been standardized initially by scaling to zero mean and unit variance.",
        "pa" = " The data has been transformed initially by setting all non-zero OTU counts to 1 (presence/absence).",
        "chi.square" = " The data has been transformed initially by applying the chi-squared transformation (Legendre & Gallagher, 2001).",
        "hellinger" = " The data has been transformed initially by applying the hellinger transformation (Legendre & Gallagher, 2001).",
        "log" = " The data has been transformed initially by applying the logarithmic transformation with logbase 2 as suggested by Anderson et al, 2006.",
        "sqrt" = " The data has been square-root transformed initially."
      )
    } else {
      " No initial data transformation has been applied."
    },
    if (any(type %in% c("pca", "ca", "pcoa", "mmds"))) {
      " The relative contribution (eigenvalue) of each axis to the total inertia in the data is indicated in percent at the axis titles."
    } else if (any(type %in% c("rda", "cca"))) {
      " The relative contribution (eigenvalue) of each axis to the total inertia in the data as well as to the constrained space only, respectively, are indicated in percent at the axis titles."
    }
  )

  # References
  references <- c(
    "1" = "Legendre, P. and Gallagher, E.D. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia 129, 271–280.",
    "2" = "Anderson, M.J., Ellingsen, K.E., McArdle, B.H. (2006). Multivariate dispersion as a measure of beta diversity. Ecology Letters 9, 683–693.",
    "3" = "Oksanen, J. (1983). Ordination of boreal heath-like vegetation with principal component analysis, correspondence analysis and multidimensional scaling. Vegetation 52, 181–189.",
    "4" = "Legendre, Pierre and Legendre, Louis (2012). Numerical Ecology. Elsevier Science. ISBN: 9780444538680.",
    "5" = "Chao, A., Chazdon, R. L., Colwell, R. K., Shen, T. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8, 148–159.",
    "6" = "Cao, Y., Williams, W.P., Bark, A.W. (1997). Effects of sample size (replicate number) on similarity measures in river benthic Aufwuchs community analysis. Water Environment Research 69, 95–106.",
    "7" = "Gower, J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics 27, 623–637.",
    "8" = "Lance, G. N. and Williams, W. T. (1966). Computer Programs for Hierarchical Polythetic Classification (“Similarity Analyses”), The Computer Journal, Volume 9, Issue 1, Pages 60–64.",
    "9" = "Bray, J. R., and Curtis, J. T. (1957). An Ordination of the upland forest community of southern Wisconsin. Ecological Monographs, 27(4), 325–349.",
    "10" = "Mountford, M. D. (1962). An index of similarity and its application to classification problems. In: P.W.Murphy (ed.), Progress in Soil Zoology, 43–50. Butterworths.",
    "11" = "Kulczynski, S. (1928). “Die Pflanzenassoziationen der Pieninen”, Bulletin international de l'Académie polonaise des Sciences et des Lettres, Classe des Sciences mathématiques et naturelles, Série B, Supplément II (1927), 57–203.",
    "12" = "Mahalanobis, PC (1936). On the generalised distance in statistics. Proceedings of the National Institute of Sciences of India 2(1): 49–55.",
    "13" = "Jaccard, P. (1901.). Étude comparative de la distribution florale dans une portion des alpes et des jura. Bulletin del la Société Vaudoise des Sciences Naturelles, 37:547-579.",
    "14" = "Morisita, M. (1959). Measuring of the dispersion and analysis of distribution patterns. Memories of the Faculty of Science, Kyushu University. Series E: Biology, 2, 215-235.",
    "15" = "Horn, H. (1966). Measurement of \"Overlap\" in Comparative Ecological Studies. The American Naturalist, 100(914), 419-424.",
    "16" = "Raup, D. M. and R. E. Crick. (1979). Measurement of faunal similarity in paleontology. Journal of Paleontology 53:1213–1227.",
    "17" = "Anderson, M.J. and Millar, R.B. (2004). Spatial variation and effects of habitat on temperate reef fish assemblages in northeastern New Zealand. Journal of Experimental Marine Biology and Ecology 305, 191–221.",
    "18" = "Clark, P. J. and Evans, F. C. (1954), Distance to Nearest Neighbor as a Measure of Spatial Relationships in Populations. Ecology, 35: 445-453. doi:10.2307/1931034."
  )

  refs2print <- integer()
  if (any(type %in% c("nmds", "mmds", "pcoa", "dca"))) {
    refs2print %<>% append(switch(distmeasure,
      "manhattan" = NULL,
      "euclidean" = NULL,
      "canberra" = 8L,
      "clark" = 18L,
      "bray" = 9L,
      "kulczynski" = 11L,
      "jaccard" = 13L,
      "gower" = 7L,
      "altGower" = 2L,
      "morisita" = 14L,
      "horn" = 15L,
      "mountford" = 10L,
      "raup" = 16L,
      "binomial" = 17L,
      "chao" = 5L,
      "cao" = 6L,
      "mahalanobis" = 12L,
      "none" = NULL
    ))
  }
  if (any(transform %in% c("hellinger", "chi.square", "log", "freq", "frequency"))) {
    refs2print %<>% append(switch(transform,
      "hellinger" = 1L,
      "chi.square" = 1L,
      "log" = 2L,
      "freq" = 3L,
      "frequency" = 3L
    ))
  }
  # combine caption and references
  captionwithrefs <- caption %>% {
    if (length(refs2print) > 0) {
      paste0(
        .,
        "\n\n\nReferences:\n\n",
        paste0(references[unique(refs2print)], collapse = "\n\n")
      )
    } else {
      return(.)
    }
  }

  # Print the caption with references, if any
  if (isTRUE(print_caption)) {
    cli::cat_line(cli::rule("Auto-generated figure caption (start)"))
    captionwithrefs %>%
      strwrap() %>%
      crayon::italic() %>%
      cli::cat_line()
    cli::cat_line(cli::rule("Auto-generated figure caption (end)"))
  }

  ##### Return  #####
  # return plot or additional details
  if (!is.null(sample_plotly)) {
    return(plotly::ggplotly(plot, tooltip = "text"))
  }
  else if (species_plotly == T) {
    return(plotly::ggplotly(plot, tooltip = "text"))
  }
  else if (!detailed_output) {
    return(plot)
  }
  else if (detailed_output) {
    if (type == "nmds") {
      screeplot <- NULL
    } else {
      ##### Screeplot  #####
      # the data for it
      if (type == "mmds" | type == "pcoa") {
        if (length(model$values$Relative_eig) > 10) {
          unconstrained_eig <- model$values$Relative_eig[1:10] * 100
        } else {
          unconstrained_eig <- model$values$Relative_eig * 100
        }
        # the scree plot
        screeplot <- ggplot(data.frame(axis = factor(as.character(c(1:length(unconstrained_eig))), levels = c(1:length(unconstrained_eig))), eigenvalues = unconstrained_eig), aes(x = axis, y = eigenvalues)) +
          geom_col() +
          # geom_text(label = round(eigenvalues, 2), vjust = -1, size = 3)  + #Can't get it to work
          theme_minimal() +
          xlab("Axis (max. 10 axes will be shown)") +
          ylab("Eigenvalue in percent of total inertia")
      } else {
        unconstrained_eig <- model$CA$eig / model$tot.chi * 100
        constrained_eig <- model$CCA$eig / model$tot.chi * 100
        if (length(constrained_eig) > 10) {
          constrained_eig <- constrained_eig[1:10]
        }
        if (length(unconstrained_eig) > 10) {
          unconstrained_eig <- unconstrained_eig[1:10]
        }
        eigenvalues <- c(constrained_eig, unconstrained_eig) # constrained combined with unconstrained
        # the scree plot
        screeplot <- ggplot(data.frame(axis = factor(names(eigenvalues), levels = names(eigenvalues)), eigenvalues = eigenvalues), aes(x = axis, y = eigenvalues)) +
          geom_col() +
          geom_text(label = round(eigenvalues, 2), vjust = -1, size = 3) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          xlab("Axis (max. 10 axes will be shown)") +
          ylab("Eigenvalue in percent of total inertia")
      }
    }

    return(list(
      plot = plot,
      screeplot = screeplot,
      figcaption = captionwithrefs,
      model = model,
      dsites = dsites,
      dspecies = dspecies,
      evf_factor_model = evf_factor_model,
      evf_numeric_model = evf_numeric_model
    ))
  }
}
