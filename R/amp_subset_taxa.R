#' Subset ampvis2 objects based on taxonomy
#'
#' Subsets the data in ampvis2 objects based on taxonomy and returns the subsetted object. See examples.
#'
#' @usage amp_subset_taxa(data, tax_vector)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param tax_vector (required) A vector with taxonomic groups, e.g. \code{c("p__Chloroflexi","p__Actinobacteria")}.
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @import dplyr
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_taxa <- function(data, tax_vector=c("p__Chloroflexi","p__Actinobacteria")) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
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
