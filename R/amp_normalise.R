#' Normalise data
#'
#' Normalises the read abundances to the total OTU counts (calculate percentages).
#'
#' @usage amp_normalise(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @import dplyr
#' @export
#' 
#' @author Rasmus Kirkegaard


amp_normalise <- function(data) {
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  # Normalise data to percent
  data$abund <- apply(data$abund,2, function(x) 100*x/sum(x)) %>% as.data.frame() 

  return(data)
}
