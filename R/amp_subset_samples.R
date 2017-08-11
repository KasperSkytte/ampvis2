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
  
  if (minreads > max(colSums(data$abund))) {
    stop(paste("Cannot subset samples with minimum", minreads, "total reads, when highest number of reads in any sample is", max(colSums(data$abund))))
  }
  
  #For printing number of removed samples
  nsamplesbefore <- nrow(d$metadata) %>% as.numeric()
  
  #Subset metadata based on ...
  data$metadata <- subset(data$metadata, ...)
  #And only keep columns in otutable that match the rows in the subsetted metadata
  data$abund <- data$abund[, rownames(data$metadata), drop = FALSE]
  
  #After subsetting the samples, remove OTU's that could possibly have 0 reads in all samples, and remove samples below the minreads
  data$abund <- data$abund[rowSums(data$abund) > 0, colSums(data$abund) >= minreads, drop = FALSE]
  #Subset the metadata again to match any removed sample(s)
  data$metadata <- data$metadata[colnames(data$abund), , drop = FALSE]
  #Subset taxonomy, and refseq if any, based on abund
  data$tax <- data$tax[rownames(data$abund), ]
  if(any(names(data) == "refseq")){
    data$refseq <- data$refseq[names(data$refseq) %in% rownames(data$abund), ]
  }
  
  nsamplesafter <- as.numeric(nrow(data$metadata))
  if (nsamplesbefore == nsamplesafter) {
    print("0 samples have been filtered.")
  } else {
    cat(paste(nsamplesbefore-nsamplesafter, "samples have been filtered \nBefore:", nsamplesbefore, "samples\nAfter:", nsamplesafter, "samples"))
  }
  
  return(data)
}
