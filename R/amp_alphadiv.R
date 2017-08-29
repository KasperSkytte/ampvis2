#' Alpha-diversity analysis
#'
#' Calculates alpha-diversity statistics for each sample and combines with the metadata.
#'
#' @usage amp_alphadiv(data, measure, rarefy)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param measure Alpha-diversity measure(s) to be included if not all. A vector of one or more of: 
#' \itemize{
#'   \item \code{"observed"}
#'   \item \code{"shannon"}
#'   \item \code{"simpson"}
#'   \item \code{"invsimpson"}
#' }
#' @param rarefy Rarefy species richness to this value. Passed directly as the \code{sample} argument to \code{\link[vegan]{rrarefy}}. (\emph{default:} \code{NULL})
#' 
#' @export
#' @import dplyr
#' @import vegan
#' 
#' @details The alpha-diversity indices are calculated per sample using the vegan function \code{\link[vegan]{diversity}}, where species richness is first rarefied using \code{\link[vegan]{rrarefy}} by the size of the \code{rarefy} argument. Refer to the vegan documentation for details about the different indices and how they are calculated. If no measure(s) are chosen, all diversity indices will be returned.
#' 
#' @examples 
#' data("AalborgWWTPs")
#' AalborgWWTPs
#' alphadiversityresult <- amp_alphadiv(AalborgWWTPs, rarefy = 20000)
#' #View(alphadiversityresult)
#' 
#' @return A data frame.
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_alphadiv <- function (data, measure = NULL, rarefy = NULL) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  abund <- data[["abund"]] %>% as.data.frame()
  Reads <- colSums(abund)
  metadata <- data[["metadata"]]
  
  if(!is.null(rarefy)){
    abund <- suppressWarnings(rrarefy(t(abund), sample = rarefy)) %>% t() %>% as.data.frame()
    if(rarefy > max(colSums(data$abund)) | rarefy < min(colSums(data$abund))) {
      warning("The chosen rarefy size is not within the min/max range of the data:\nMin#reads: ", as.character(min(colSums(data$abund))), "\nMax#reads: ", as.character(max(colSums(data$abund))))
    }
    if (min(colSums(abund)) < rarefy) {
      warning("The following samples have not been rarefied (less than ", as.character(rarefy), " reads):\n", paste(rownames(metadata[which(colSums(abund) < rarefy),]), collapse = ", "))
    }
  }
  
  if (!any(abund == 1)) {
    warning("The data you have provided does not have\n", 
            "any singletons. This is highly suspicious. Results of richness\n", 
            "estimates (for example) are probably unreliable, or wrong, if you have already\n", 
            "trimmed low-abundance taxa from the data.\n", "\n", 
            "We recommend that you find the un-trimmed data and retry.")
  }
  
  abund <- as.matrix(t(abund))
  
  renamevec = c("observed", "chao1", "ace", "shannon", "simpson", 
                "invsimpson")
  
  names(renamevec) <- c("s.obs", "s.chao1", "s.ace", "shannon", 
                        "simpson", "invsimpson")
  
  if (is.null(measure)) {
    measure = renamevec %>% as.character() %>% tolower()
  }
  if (any(measure %in% names(renamevec))) {
    measure[measure %in% names(renamevec)] <- renamevec[names(renamevec) %in% measure]
  }
  if (!any(measure %in% renamevec)) {
    stop("None of the `measure`s provided are supported. Try default `NULL` instead.")
  }
  
  outlist = vector("list")
  
  estimRmeas = c("chao1", "observed", "ace")
  
  if (any(estimRmeas %in% measure)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(abund)))))
  }
  if ("shannon" %in% measure) {
    outlist <- c(outlist, list(shannon = vegan::diversity(abund, index = "shannon")))
  }
  if ("simpson" %in% measure) {
    outlist <- c(outlist, list(simpson = vegan::diversity(abund, index = "simpson")))
  }
  if ("invsimpson" %in% measure) {
    outlist <- c(outlist, list(invsimpson = vegan::diversity(abund, index = "invsimpson")))
  }
  
  rich = do.call("cbind", outlist)
  namechange = intersect(colnames(rich), names(renamevec))
  
  colnames(rich)[colnames(rich) %in% namechange] <- renamevec[namechange]
  
  colkeep = sapply(paste0("(se\\.){0,}", measure), grep, 
                   colnames(rich), ignore.case = TRUE)
  
  rich = rich[, sort(unique(unlist(colkeep))), drop = FALSE]
  
  #Combine metadata and results and return
  combined <- data.frame(metadata, Reads = Reads, rich) %>% arrange(Reads)
  return(combined)
}

