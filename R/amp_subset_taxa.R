#' Subset data based on taxonomy.
#'
#' Subset data based on taxonomy.
#'
#' @usage amp_subset_taxa(data, ...)
#'
#' @param data (required) A object.
#' @param ... (required) Additional data is passed on to the R subset function. The samples can be subset based on any variable in the taxonomy.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_taxa <- function(data, ...) {
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
  
  #extract data from the list
  metadata <- data$metadata
  abund <- data$abund
  tax <- data$tax
  

  #subset tax table based on ... and only keep rows in abund and metadata matching the rows in the subsetted tax table
  #data$tax <- subset(tax, ...)
  
  data$tax <- subset(tax, Kingdom == "k__Bacteria")
  
  data$abund <- abund[as.character(data$tax$OTU),]
  #data$abund <- subset(data$abund, data$abund$OTU %in% data$tax$OTU)
  #data$metadata <- metadata[colnames(data$abund), , drop=FALSE]
  
  if(any(names(data) == "refseq")){
    data$refseq <- data$refseq[names(data$refseq) %in% rownames(data$abund), ]
  }
  
  return(data)
}
