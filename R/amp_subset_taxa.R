#' Subset ampvis2 objects based on taxonomy
#'
#' Subsets the data in ampvis2 objects based on taxonomy and returns the subsetted object. 
#'
#' @usage amp_subset_taxa(data, tax_vector)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_vector (required) A vector with the taxonomic groups with which to perform the subset. The prefixes \code{"k__"}, \code{"p__"}, \code{"c__"}, \code{"o__"}, \code{"f__"}, \code{"g__"}, \code{"s__"} indicates their taxonomic rank and the following characters their name (almost always with a capital first letter), e.g. \code{c("p__Chloroflexi","p__Actinobacteria")}.
#' @param normalise (\emph{logical}) Normalise the read abundances to the total amount of reads (percentages) \emph{BEFORE} the subset. (\emph{default:} \code{FALSE})
#' @param remove (\emph{logical}) If set to TRUE, then the taxa matching the provided vector will be removed instead of being the only ones kept in the data. (\emph{default:} \code{FALSE})
#' 
#' @return A modifed ampvis2 object
#' @import ape
#' @export
#' 
#' @details 
#' The taxonomy subset is done by providing a \code{tax_vector} of taxa names which are then matched to the taxonomy table, where all other taxa not matching the \code{tax_vector} are removed. If \code{remove = TRUE}, then the matching taxa are the ones being removed instead. The taxa names in \code{tax_vector} will be matched in all columns of the taxonomy table.
#' 
#' @section Normalising data for use in heatmaps:
#' By default the raw read counts in the abundance matrix are normalised (transformed to percentages) by \code{\link{amp_heatmap}} automatically. This means that the relative abundances shown will be calculated based on the remaining taxa after the subset, not including the removed taxa, if any. To circumvent this, set \code{normalise = TRUE} when subsetting with the \code{\link{amp_subset_taxa}} and \code{\link{amp_subset_samples}} functions and then set \code{raw = TRUE} when using \code{\link{amp_heatmap}}, see the example below.
#' 
#' \preformatted{
#' data("AalborgWWTPs")
#' subsettedData <- amp_subset_taxa(AalborgWWTPs,
#'                                  tax_vector = c("p__Chloroflexi", "p__Actinobacteria"),
#'                                  normalise = TRUE
#'                                  )
#' amp_heatmap(subsettedData,
#'             group_by = "Plant",
#'             tax_aggregate = "Phylum",
#'             tax_add = "Genus",
#'             raw = TRUE
#'             )
#' }
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #An overview heatmap of the data:
#' amp_heatmap(AalborgWWTPs, 
#'             tax_aggregate = "Genus",
#'             group_by = "Plant",
#'             tax_add = "Phylum"
#'             )
#' 
#' #Remove all taxa except the phyla Chloroflexi and Actinobacteria 
#' #and the Genera Rhodoferax and Trichococcus:
#' tax_vector <- c("p__Chloroflexi",
#'             "p__Actinobacteria",
#'             "g__Rhodoferax",
#'             "g__Trichococcus")
#'             
#' AalborgWWTPs_subset <- amp_subset_taxa(AalborgWWTPs, 
#'                                        tax_vector = tax_vector
#'                                        )
#' 
#' #The resulting subset:
#' amp_heatmap(AalborgWWTPs_subset, 
#'             tax_aggregate = "Genus",
#'             group_by = "Plant",
#'             tax_add = "Phylum"
#'             )
#' 
#' #Or if remove = TRUE then the taxa in tax_vector are the ones being removed:
#' AalborgWWTPs_subset <- amp_subset_taxa(AalborgWWTPs, 
#'                                        tax_vector = tax_vector,
#'                                        remove = TRUE
#'                                        )
#' #The resulting subset:
#' amp_heatmap(AalborgWWTPs_subset, 
#'             tax_aggregate = "Genus",
#'             group_by = "Plant",
#'             tax_add = "Phylum"
#'             )
#' 
#' @seealso 
#' \code{\link{amp_subset_samples}}, \code{\link{amp_heatmap}}
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @author Rasmus Hansen Kirkegaard \email{rhk@@bio.aau.dk}


amp_subset_taxa <- function(data, 
                            tax_vector = c("p__Chloroflexi", "p__Actinobacteria"), 
                            normalise = FALSE,
                            remove = FALSE) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ### Check if refseq data is in the right format
  if(!is.null(data$refseq) & !class(data$refseq) == "DNAbin") {
    stop("The refseq element is not of class \"DNAbin\". The reference sequences must be loaded with ape::read.dna().")
  }
  
  ### calculate percentages 
  if (normalise) {
    data$abund <- as.data.frame(apply(data$abund,2, function(x) 100*x/sum(x)))
  }
  
  # Make new list
  selection <- c(which(data$tax$Kingdom %in% tax_vector),
                 which(data$tax$Phylum %in% tax_vector),
                 which(data$tax$Class %in% tax_vector),
                 which(data$tax$Order %in% tax_vector),
                 which(data$tax$Family %in% tax_vector),
                 which(data$tax$Genus %in% tax_vector),
                 which(data$tax$Species %in% tax_vector),
                 which(data$tax$OTU %in% tax_vector)
  )
  selection <- unique(selection)
  newtax <- data$tax[selection,]
  if(isTRUE(remove)) {
    data$tax <- subset(data$tax, !OTU %in% newtax$OTU)
  } else if(!isTRUE(remove)) {
    data$tax <- newtax
  }
  data$abund <- data$abund[rownames(data$tax),]
  
  if (any(names(data) == "refseq")) {
    data$refseq <- data$refseq[rownames(data$tax)]
  }
  return(data)
}
