#' Subset data based on taxonomy.
#'
#' Subset data based on taxonomy.
#'
#' @usage amp_subset_taxa(data, tax_vector)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param tax_vector (required) A vector with taxonomic groups, e.g. \code{c("p__Chloroflexi","p__Actinobacteria")}
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_taxa <- function(data, tax_vector=c("p__Chloroflexi","p__Actinobacteria")) {
  #Check the data first
  if(!is.list(data) | 
     !any(names(data) == "abund") |
     !any(names(data) == "tax") | 
     !any(names(data) == "metadata") | 
     !is.data.frame(data[["abund"]]) |
     !is.data.frame(data[["tax"]]) |
     !is.data.frame(data[["metadata"]])
  ) {
    stop("The data must be a list with three dataframes named abund, tax and metadata")
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
