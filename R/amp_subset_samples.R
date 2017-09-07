#' Subset ampvis2 objects based on sample metadata
#'
#' Subsets the data in ampvis2 objects based on metadata and returns the subsetted object. 
#'
#' @usage amp_subset_samples(data, ...)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param minreads Minimum number of reads pr. sample. (\emph{default:} \code{1})
#' @param ... Logical expression indicating elements or rows to keep in the metadata. Missing values are taken as false. Directly passed to \code{subset()}. 
#' @param normalise (\emph{logical}) Normalise the read abundances to the total amount of reads (percentages) \emph{BEFORE} the subset. (\emph{default:} \code{FALSE})
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @import dplyr
#' @import ape
#' @import stringr
#' @export
#' 
#' @details The subset is performed on the metadata by \code{subset()} and the abundance- and taxonomy tables are then adjusted accordingly.
#' 
#' @examples 
#' #Load example data
#' data("MiDAS")
#' 
#' #Show a short summary about the data by simply typing the name of the object in the console
#' MiDAS
#' 
#' #Keep only samples containing Aalborg West or East in the Plant column
#' MiDASsubset <- amp_subset_samples(MiDAS, Plant %in% c("Aalborg West", "Aalborg East"))
#' 
#' #Summary
#' MiDASsubset
#' 
#' #Keep only samples containing Aalborg West or East in the Plant column 
#' #and remove the sample "16SAMP-749". Remove any sample(s) with less than 20000 total reads
#' MiDASsubset2 <- amp_subset_samples(MiDAS,
#'     Plant %in% c("Aalborg West", "Aalborg East") & !SampleID %in% c("16SAMP-749"), 
#'     minreads = 20000)
#'     
#' #Summary
#' MiDASsubset2
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_samples <- function(data, ..., minreads = 1, normalise = FALSE) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  if (minreads > max(colSums(data$abund))) {
    stop(paste("Cannot subset samples with minimum", minreads, "total reads, when highest number of reads in any sample is", max(colSums(data$abund))))
  }
  
  ### Check if refseq data is in the right format
  if(!is.null(data$refseq) & !class(data$refseq) == "DNAbin") {
    stop("The refseq element is not of class \"DNAbin\". The reference sequences must be loaded with ape::read.dna().")
  }
  
  ### calculate percentages 
  if (normalise) {
    data$abund <- apply(data$abund,2, function(x) 100*x/sum(x)) %>% as.data.frame() 
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
  
  #Subset taxonomy based on abund
  data$tax <- data$tax[rownames(data$abund), ]
  
  #Subset refseq, if any, based on abund
  if(any(names(data) == "refseq")){
    #sometimes there is taxonomy alongside the OTU ID's. Anything after a ";" will be ignored
    names_stripped <- stringr::str_split(names(data$refseq), ";", simplify = TRUE)[,1]
    data$refseq <- data$refseq[names_stripped %in% rownames(data$abund)] 
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
