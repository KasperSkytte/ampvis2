#' Subset data based on metadata.
#'
#' Subset data based on metadata.
#'
#' @usage amp_subset_samples(data, ...)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param minreads (optional) Minimum number of reads pr. sample. (\emph{default:} \code{1})
#' @param ... (optional) Additional arguments for the \code{subset()} function. 
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
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
  
  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(data$metadata) %>% as.numeric()
  nOTUsbefore <- nrow(data$abund) %>% as.numeric()
  
  #Subset metadata based on ...
  data$metadata <- subset(data$metadata, ...)
  data$metadata <- droplevels(data$metadata) #Drop unused factor levels or fx heatmaps will show a "NA" column
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
  
  #Print number of removed samples
  nsamplesafter <- nrow(data$metadata) %>% as.numeric()
  nOTUsafter <- nrow(data$abund) %>% as.numeric()
  if (nsamplesbefore == nsamplesafter) {
    print("0 samples have been filtered.")
  } else {
    cat(paste(nsamplesbefore-nsamplesafter, "samples and", nOTUsbefore-nOTUsafter,"OTUs have been filtered \nBefore:", nsamplesbefore, "samples and", nOTUsbefore, "OTUs\nAfter:", nsamplesafter, "samples and", nOTUsafter, "OTUs"))
  }
  
  return(data)
}
