#' Subset ampvis2 objects based on taxonomy
#'
#' Subsets the data in ampvis2 objects based on taxonomy and returns the subsetted object. See examples.
#'
#' @usage amp_subset_taxa(data, tax_vector)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_vector (required) A vector with taxonomic groups, e.g. \code{c("p__Chloroflexi","p__Actinobacteria")}.
#' @param normalise (\emph{logical}) Normalise the read abundances to the total amount of reads (percentages) \emph{BEFORE} the subset. (\emph{default:} \code{FALSE})
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @import ape
#' @export
#' 
#' @details 
#' The subset of taxa is done by providing a vector of taxa names, where all other taxa are then removed. The taxa names in the vector will be matched in all columns in the taxonomy table.
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
#' data("MiDAS")
#' 
#' #Show a short summary about the data by simply typing the name of the object in the console
#' MiDAS
#' 
#' #Remove all taxa except the phyla Chloroflexi and Actinobacteria
#' MiDASsubset <- amp_subset_taxa(MiDAS, tax_vector = c("p__Chloroflexi", "p__Actinobacteria"))
#' 
#' #Summary
#' MiDASsubset
#' 
#' @seealso 
#' \code{\link{amp_subset_samples}}, \code{\link{amp_heatmap}}
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @author Rasmus Hansen Kirkegaard \email{rhk@bio.aau.dk}


amp_subset_taxa <- function(data, tax_vector = c("p__Chloroflexi", "p__Actinobacteria"), normalise = FALSE) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ### Check if refseq data is in the right format
  if(!is.null(data$refseq) & !class(data$refseq) == "DNAbin") {
    stop("The refseq element is not of class \"DNAbin\". The reference sequences must be loaded with ape::read.dna().")
  }
  
  ### calculate percentages 
  if (normalise) {
    data$abund <- apply(data$abund,2, function(x) 100*x/sum(x)) %>% as.data.frame() 
  }
  
  # Make selection
  selection <- c(which(data$tax$Kingdom %in% tax_vector),
                 which(data$tax$Phylum %in% tax_vector),
                 which(data$tax$Class %in% tax_vector),
                 which(data$tax$Order %in% tax_vector),
                 which(data$tax$Family %in% tax_vector),
                 which(data$tax$Genus %in% tax_vector),
                 which(data$tax$Species %in% tax_vector),
                 which(data$tax$OTU %in% tax_vector)
                 )
  
  # Make new list
  data$abund <- data$abund[selection,]
  data$tax <- data$tax[selection,]
  data$metadata <- data$metadata
  
  if (any(names(data) == "refseq")) {
    data$refseq <- data$refseq[selection]
  }

  return(data)
}
