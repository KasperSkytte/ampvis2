#' Subset data based on taxonomy.
#'
#' Subset data based on taxonomy.
#'
#' @usage amp_subset_taxa(data, ...)
#'
#' @param data (required) A object.
#' @param tax_vector (required) a vector with taxonomic groups e.g. c("p__Chloroflexi","p__Actinobacteria")
#' 
#' @return A list object.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_taxa<- function(data, tax_vector=c("p__Chloroflexi","p__Actinobacteria")) {
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
  #selection<-which(data$tax$...)
  selection<-c(which(data$tax$Kingdom %in% tax_vector),
               which(data$tax$Phylum %in% tax_vector),
               which(data$tax$Class %in% tax_vector),
               which(data$tax$Order %in% tax_vector),
               which(data$tax$Family %in% tax_vector),
               which(data$tax$Species %in% tax_vector),
               which(data$tax$OTU %in% tax_vector)
               )
  
  #newDF <- subset(oldDF, ...)
  
  # Make new list
  d1<-data$abund[selection,]
  d2<-data$tax[selection,]
  d3<-data$metadata
  d4<-data$refseq[selection]
  data<- list(abund = d1, tax = d2, metadata = d3, refseq = d4)

  return(data)
}
