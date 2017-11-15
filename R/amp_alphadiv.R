#' Alpha-diversity analysis
#'
#' Calculate alpha-diversity indices for each sample and combines with the metadata.
#'
#' @usage amp_alphadiv(data, measure, rarefy)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param measure Alpha-diversity measure(s) to be included if not all. A vector of one or more of: 
#' \itemize{
#'   \item \code{"observed"}
#'   \item \code{"shannon"}
#'   \item \code{"simpson"}
#'   \item \code{"invsimpson"}
#' }
#' @param richness (\emph{logical}) Also calculate sample richness estimates (Chao1 and ACE) as calculated by \code{\link[vegan]{estimateR}}. (\emph{default:} \code{FALSE})
#' @param rarefy Rarefy species richness to this value before calculating alpha diversity and/or richness. Passed directly as the \code{sample} argument to \code{\link[vegan]{rrarefy}}. (\emph{default:} \code{NULL})
#' 
#' @export
#' @import dplyr
#' @import vegan
#' 
#' @details The alpha-diversity indices are calculated per sample using the vegan function \code{\link[vegan]{diversity}}, where the read abundances are first rarefied using \code{\link[vegan]{rrarefy}} by the size of the \code{rarefy} argument. Refer to the vegan documentation for details about the different indices and how they are calculated. If no measure(s) are chosen, all diversity indices will be returned.
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Subsample/rarefy to 20000 reads and then calculate Shannon and Simpson alpha-diversity indices
#' alphadiversityresult <- amp_alphadiv(AalborgWWTPs, measure = c("shannon", "simpson"), rarefy = 20000)
#' 
#' #Explore the results in the data frame
#' #View(alphadiversityresult)
#' 
#' @return A data frame.
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_alphadiv <- function (data,
                          measure = NULL,
                          richness = FALSE,
                          rarefy = NULL) {
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  #check measures
  validMeasures <- c("observed", "shannon", "simpson", "invsimpson")
  if(is.null(measure)) {
    measure <- validMeasures
  } else if(!is.null(measure) & any(!measure %in% validMeasures)) {
    measure <- measure %>% tolower()
    warning("Some or none of the provided measures were not recognised, calculating all. Valid options are:\n", paste0(validMeasures, collapse = ", "))
    measure <- validMeasures
  }
  
  results <- data$metadata
  names <- results[,1] %>% as.vector() #for making sure the ordering of the values in the vectors calculated later match the order of the metadata samples
  abund <- t(data$abund)
  
  #Add Reads column
  Reads <- colSums(data$abund)
  Reads <- Reads[names]
  results$Reads <- Reads
  
  if(!is.null(rarefy) & is.numeric(rarefy)){
    if(rarefy > max(results$Reads) ) {
      stop("The chosen rarefy size is larger than the largest amount of reads in any sample (", as.character(max(results$Reads)), ").")
    } else if (rarefy < min(results$Reads)) {
      abund <- suppressWarnings(vegan::rrarefy(abund, sample = rarefy)) %>% as.data.frame()
      warning("The chosen rarefy size (", as.character(rarefy), ") is smaller than the smallest amount of reads in any sample (", as.character(min(colSums(data$abund))), ").")
    } else {
      abund <- suppressWarnings(vegan::rrarefy(abund, sample = rarefy)) %>% as.data.frame()
      if (min(results$Reads) < rarefy) {
        message("The following samples have not been rarefied (less than ", as.character(rarefy), " reads):\n", paste(rownames(data$metadata[which(results$Reads < rarefy),]), collapse = ", "))
      }
    }
  } else if(!is.null(rarefy) & !is.numeric(rarefy)) {
    stop("Argument rarefy must be numerical.")
  }
  
  #warning from phyloseq::estimate_richness
  if (!any(abund == 1)) {
    warning("The data you have provided does not have\n", 
            "any singletons. This is highly suspicious. Results of richness\n", 
            "estimates (for example) are probably unreliable, or wrong, if you have already\n", 
            "trimmed low-abundance taxa from the data.\n", "\n", 
            "We recommend that you find the un-trimmed data and retry.")
  }
  
  if(any("observed" %in% measure) | is.null(measure)) {
    ObservedOTUs <- colSums(t(abund) > 0)
    results$ObservedOTUs <- ObservedOTUs[names]
  }
  if(any("shannon" %in% measure)) {
    Shannon <- vegan::diversity(abund, index = "shannon")
    results$Shannon = Shannon[names]
  }
  if(any("simpson" %in% measure)) {
    Simpson <- vegan::diversity(abund, index = "simpson")
    results$Simpson <- Simpson[names]
  }
  if(any("invsimpson" %in% measure)) {
    invSimpson <- vegan::diversity(abund, index = "invsimpson")
    results$invSimpson <- invSimpson[names]
  }
  if(richness) {
    richness <- t(vegan::estimateR(abund)) %>% as.data.frame()
    richness <- richness[names,, drop = FALSE]
    results$Chao1 <- richness[,"S.chao1"]
    results$ACE <- richness[,"S.ACE"]
  }
  
  results <- results %>% arrange(Reads)
  return(results)
}
