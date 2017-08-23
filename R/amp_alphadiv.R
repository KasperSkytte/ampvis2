#' Alpha-diversity analysis
#'
#' Calculates alpha-diversity statistics for each sample.
#'
#' @usage amp_alphadiv(data, measure = c(""), rarefy = 10000)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param measure Alpha-diversity measure(s) to be included if not all. A vector of one or more of: 
#' \itemize{
#'   \item \code{"observed"}
#'   \item \code{"shannon"}
#'   \item \code{"simpson"}
#'   \item \code{"invsimpson"}
#' }
#' @param rarefy Rarefy species richness to this value. (\emph{default:} \code{10000})
#' 
#' @export
#' @import tidyverse
#' @import vegan
#' @details \code{amp_alphadiv} simply calculates alpha-diversity indices using the vegan function \code{\link[vegan]{diversity}} and returns the results combined with the metadata from the provided ampvis2 object. For more details about the exact formulas used see the references.
#' @references 
#'   \url{https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf}
#'   
#'   Hill, M. (1973). "Diversity and Evenness: A Unifying Notation and Its Consequences". Ecology, 54(2), 427-432. \url{doi.org/10.2307/1934352}
#' @examples
#' 
#' @return A data frame.
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_alphadiv <- function (data, measure = NULL, rarefy = 10000) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  abund <- data[["abund"]] %>% as.data.frame()
  Reads <- colSums(abund)
  metadata <- data[["metadata"]]

  
  if(!is.null(rarefy)){
    abund <- rrarefy(t(abund), sample = rarefy) %>% t() %>% as.data.frame()
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
  
  combined <- data.frame(metadata, Reads = Reads, rich) %>% arrange(Reads)
  return(combined)
}

