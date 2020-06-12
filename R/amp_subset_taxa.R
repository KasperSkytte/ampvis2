#' Subset ampvis2 objects based on taxonomy
#'
#' Subsets the data in ampvis2 objects based on taxonomy and returns the subsetted object.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_vector A character vector with the exact names of taxa to keep. This vector is matched as-is on all taxonomic ranks, so remember to use prefixes if used in your taxonomy, e.g. \code{c("p__Chloroflexi","p__Actinobacteria")}. (\emph{default:} \code{NULL})
#' @param normalise (\emph{logical}) Normalise the OTU read counts to 100 (ie percent) per sample \emph{BEFORE} the subset. (\emph{default:} \code{FALSE})
#' @param remove (\emph{logical}) If set to TRUE, then the taxa matching the provided vector will be removed instead of being the only ones kept in the data. (\emph{default:} \code{FALSE})
#'
#' @return A modifed ampvis2 object
#' @importFrom ape drop.tip
#' @export
#'
#' @details
#' The taxonomy subset is done by providing a \code{tax_vector} of taxa names which are then matched to the taxonomy table, where all other taxa not matching the \code{tax_vector} are removed. If \code{remove = TRUE}, then the matching taxa are the ones being removed instead. The taxa names in \code{tax_vector} will be matched in all columns of the taxonomy table.
#'
#' @section Preserving relative abundances in a subset of larger data:
#' By default the raw read counts in the abundance matrix are normalised (transformed to percentages) by some plotting functions automatically (for example \code{\link{amp_heatmap}}, \code{\link{amp_timeseries}}, and more). This means that the relative abundances shown will be calculated based on the remaining taxa after the subset, not including the removed taxa, if any. To circumvent this, set \code{normalise = TRUE} when subsetting with the \code{\link{amp_subset_taxa}} and \code{\link{amp_subset_samples}} functions, and then set \code{raw = TRUE} in the plotting function. This will transform the OTU counts to relative abundances BEFORE the subset, and setting \code{raw = TRUE} will skip the transformation in the plotting function, see the example below.
#'
#' \preformatted{
#' data("MiDAS")
#' subsettedData <- amp_subset_samples(MiDAS,
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
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # An overview heatmap of the data:
#' amp_heatmap(AalborgWWTPs,
#'   tax_aggregate = "Genus",
#'   group_by = "Plant",
#'   tax_add = "Phylum"
#' )
#'
#' # Remove all taxa except the phyla Chloroflexi and Actinobacteria
#' # and the Genera Rhodoferax and Trichococcus:
#' tax_vector <- c(
#'   "p__Chloroflexi",
#'   "p__Actinobacteria",
#'   "g__Rhodoferax",
#'   "g__Trichococcus"
#' )
#'
#' AalborgWWTPs_subset <- amp_subset_taxa(AalborgWWTPs,
#'   tax_vector = tax_vector
#' )
#'
#' # The resulting subset:
#' amp_heatmap(AalborgWWTPs_subset,
#'   tax_aggregate = "Genus",
#'   group_by = "Plant",
#'   tax_add = "Phylum"
#' )
#'
#' # Or if remove = TRUE then the taxa in tax_vector are the ones being removed:
#' AalborgWWTPs_subset <- amp_subset_taxa(AalborgWWTPs,
#'   tax_vector = tax_vector,
#'   remove = TRUE
#' )
#' # The resulting subset:
#' amp_heatmap(AalborgWWTPs_subset,
#'   tax_aggregate = "Genus",
#'   group_by = "Plant",
#'   tax_add = "Phylum"
#' )
#' @seealso
#' \code{\link{amp_load}}, \code{\link{amp_subset_samples}}
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @author Rasmus Hansen Kirkegaard \email{rhk@@bio.aau.dk}

amp_subset_taxa <- function(data,
                            tax_vector = NULL,
                            normalise = FALSE,
                            remove = FALSE) {

  ### Data must be in ampvis2 format
  is_ampvis2(data)

  if (!is.null(data$refseq)) {
    if (!class(data$refseq) == "DNAbin") {
      stop("The refseq element is not of class \"DNAbin\". The reference sequences must be loaded with ape::read.dna().", call. = FALSE)
    }
  }

  if (!is.null(data$tree)) {
    if (!class(data$tree) == "phylo") {
      stop("The tree element is not of class \"phylo\". The tree must be loaded with ape::read.tree().", call. = FALSE)
    }
  }

  # For printing removed OTUs
  nOTUsbefore <- nrow(data$abund)

  # normalise counts before subsetting
  if (isTRUE(normalise)) {
    data <- normaliseTo100(data)
  }
  
  if(!is.null(tax_vector)) {
    # Match tax_vector with all taxonomic levels
    selection <- c(
      which(data$tax$Kingdom %in% tax_vector),
      which(data$tax$Phylum %in% tax_vector),
      which(data$tax$Class %in% tax_vector),
      which(data$tax$Order %in% tax_vector),
      which(data$tax$Family %in% tax_vector),
      which(data$tax$Genus %in% tax_vector),
      which(data$tax$Species %in% tax_vector),
      which(data$tax$OTU %in% tax_vector)
    )
    selection <- unique(selection)
    newtax <- data$tax[selection, ]
  
    # subset
    if (isTRUE(remove)) {
      data$tax <- subset(data$tax, !OTU %in% newtax$OTU)
    } else if (!isTRUE(remove)) {
      data$tax <- newtax
    }
    data$abund <- data$abund[rownames(data$abund) %in% rownames(data$tax), , drop = FALSE]
  
    # subset sequences
    if (any(names(data) == "refseq")) {
      data$refseq <- data$refseq[rownames(data$tax)]
    }
  
    # subset tree
    if (!is.null(data$tree)) {
      data$tree <- ape::drop.tip(
        phy = data$tree,
        tip = data$tree$tip.label[!data$tree$tip.label %in% data$tax$OTU]
      )
    }
  }

  if(any(colSums(data$abund) == 0))
    warning("One or more samples have 0 total reads", call. = FALSE)
  nOTUsafter <- nrow(data$abund)
  if (nOTUsbefore == nOTUsafter) {
    message("0 OTU's have been filtered.")
  } else {
    message(paste(nOTUsbefore - nOTUsafter, "OTUs have been filtered \nBefore:", nOTUsbefore, "OTUs\nAfter:", nOTUsafter, "OTUs"))
  }

  return(data)
}
