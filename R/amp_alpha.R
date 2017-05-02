#' Calculate alpha-diversity statistics for each sample
#'
#' Calculate alpha-diversity statistics for each sample and combine it with metadata.
#'
#' @usage amp_alpha(data)
#'
#' @param data (required) An ampvis object.
#' @param measure Alpha-diversity measures to be included (default:Observed).
#' 
#' @export
#' @import dplyr
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_alpha <- function (data, measures = NULL, rarefy = 10000) {

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
            "We recommended that you find the un-trimmed data and retry.")
  }
  
  abund <- as.matrix(t(abund))
  
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
                "InvSimpson")
  
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
                        "simpson", "invsimpson")
  
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  
  outlist = vector("list")
  
  estimRmeas = c("Chao1", "Observed", "ACE")
  
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(abund)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = vegan::diversity(abund, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = vegan::diversity(abund, index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = vegan::diversity(abund, index = "invsimpson")))
  }
  
  rich = do.call("cbind", outlist)
  namechange = intersect(colnames(rich), names(renamevec))
  colnames(rich)[colnames(rich) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, 
                   colnames(rich), ignore.case = TRUE)
  rich = rich[, sort(unique(unlist(colkeep))), drop = FALSE]
  rich <- as.data.frame(rich)
  combined <- cbind.data.frame(metadata, Reads, rich) %>% arrange(Reads)
  return(combined)
}

