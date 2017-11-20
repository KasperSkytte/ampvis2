#' Subset ampvis2 objects based on sample metadata
#'
#' Subsets the data in ampvis2 objects based on metadata and returns the subsetted object. 
#'
#' @usage amp_subset_samples(data, ...)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param ... Logical expression indicating elements or rows to keep in the metadata. Missing values are taken as false. Directly passed to \code{subset()}. 
#' @param minreads Minimum number of reads pr. sample. Samples below this value will be removed. (\emph{default:} \code{0})
#' @param normalise (\emph{logical}) Normalise the read abundances to the total amount of reads (percentages) \emph{BEFORE} the subset. (\emph{default:} \code{FALSE})
#' @param removeAbsents (\emph{logical}) Whether to remove OTU's that may have 0 read abundance in all samples after the subset. (\emph{default:} \code{TRUE})
#' 
#' @return A modifed ampvis2 object
#' @import dplyr
#' @import ape
#' @import stringr
#' @export
#' 
#' @details The subset is performed on the metadata by \code{subset()} and the abundance- and taxonomy tables are then adjusted accordingly.
#' 
#' @section Normalising data for use in heatmaps:
#' By default the raw read counts in the abundance matrix are normalised (transformed to percentages) by \code{\link{amp_heatmap}} automatically. This means that the relative abundances shown will be calculated based on the remaining taxa after the subset, not including the removed taxa, if any. To circumvent this, set \code{normalise = TRUE} when subsetting with the \code{\link{amp_subset_taxa}} and \code{\link{amp_subset_samples}} functions and then set \code{raw = TRUE} when using \code{\link{amp_heatmap}}, see the example below.
#' 
#' \preformatted{
#' data("MiDAS")
#' subsettedData <- amp_subset_samples(MiDAS,
#'                                     Plant \%in\% c("Aalborg West", "Aalborg East"),
#'                                     normalise = TRUE
#'                                     )
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
#' @seealso 
#' \code{\link{amp_subset_taxa}}, \code{\link{amp_heatmap}}
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_subset_samples <- function(data, ..., minreads = 0, normalise = FALSE, removeAbsents = TRUE) {
  
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
  
  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(data$metadata) %>% as.numeric()
  nOTUsbefore <- nrow(data$abund) %>% as.numeric()
  
  #remove samples below minreads BEFORE percentages
  data$abund <- data$abund[, colSums(data$abund) >= minreads, drop = FALSE]
  
  #Subset the metadata again to match any removed sample(s)
  data$metadata <- data$metadata[which(rownames(data$metadata) %in% colnames(data$abund)), , drop = FALSE]
  
  ### calculate percentages 
  if (normalise == TRUE) {
    data$abund <- apply(data$abund,2, function(x) 100*x/sum(x)) %>% as.data.frame() 
  }
  
  #Subset metadata based on ...
  data$metadata <- subset(data$metadata, ...)
  data$metadata <- droplevels(data$metadata) #Drop unused factor levels or fx heatmaps will show a "NA" column
  
  #And only keep columns in otutable that match the rows in the subsetted metadata
  data$abund <- data$abund[, which(colnames(data$abund) %in% rownames(data$metadata)), drop = FALSE]
  
  #After subsetting the samples, remove OTU's that may have 0 reads in all samples
  if(removeAbsents == TRUE) {
    data$abund <- data$abund[rowSums(data$abund) > 0, , drop = FALSE]
  }
  
  #Subset taxonomy based on abund
  data$tax <- data$tax[which(rownames(data$tax) %in% rownames(data$abund)), , drop = FALSE]
  
  #make sure the order of sample names are identical between abund and metadata
  data$abund = data$abund[,rownames(data$metadata), drop = FALSE]
  data$tax = data$tax[rownames(data$abund),, drop = FALSE]
  
  #Subset refseq, if any, based on abund
  if(any(names(data) == "refseq")){
    if(!is.null(names(data$refseq))) {
      #sometimes there is taxonomy alongside the OTU ID's. Anything after a ";" will be ignored
      names_stripped <- stringr::str_split(names(data$refseq), ";", simplify = TRUE)[,1]
      data$refseq <- data$refseq[names_stripped %in% rownames(data$abund)] 
    } else if(is.null(names(data$refseq))) {
      warning("DNA sequences have not been subsetted, could not find the names of the sequences in data$refseq.")
    }
  }
  
  #Print number of removed samples
  nsamplesafter <- nrow(data$metadata) %>% as.numeric()
  nOTUsafter <- nrow(data$abund) %>% as.numeric()
  if (nsamplesbefore == nsamplesafter) {
    message("0 samples have been filtered.")
  } else {
    message(paste(nsamplesbefore-nsamplesafter, "samples and", nOTUsbefore-nOTUsafter,"OTUs have been filtered \nBefore:", nsamplesbefore, "samples and", nOTUsbefore, "OTUs\nAfter:", nsamplesafter, "samples and", nOTUsafter, "OTUs"))
  }
  
  return(data)
}
