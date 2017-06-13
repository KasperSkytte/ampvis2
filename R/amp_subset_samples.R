#' Subset data based on metadata.
#'
#' Subset data based on metadata.
#'
#' @usage amp_subset_samples(data, ...)
#'
#' @param data (required) A object.
#' @param minreads (optional) Minimum number of reads pr. sample (default: 1).
#' @param ... (optional) Additional data is passed on to the R subset function. The samples can be subset based on any variable in the associated metadata.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_samples <- function(data, ..., minreads = 1) {
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

  goodSamples <- data.frame(SeqID = colnames(abund), Seqs = colSums(abund)) %>% subset(Seqs >= minreads)
  
  metadata <- subset(metadata, SeqID %in% goodSamples$SeqID)
  
  #subset metadata based on ... and only keep columns in otutable matching the rows in the subsetted metadata
  data$metadata <- subset(metadata, ...)
  data$metadata <- droplevels(data$metadata)
  data$abund <- abund[, data$metadata$SeqID, drop=FALSE]
  data$abund <- data$abund[!apply(data$abund, 1, function(row) all(row <= 0)),] #remove low abundant OTU's 
  data$tax <- data$tax[rownames(data$abund),] #same with taxonomy
  
  if(any(names(data) == "refseq")){
    data$refseq <- data$refseq[names(data$refseq) %in% rownames(data$abund), ]
  }

  return(data)
}
