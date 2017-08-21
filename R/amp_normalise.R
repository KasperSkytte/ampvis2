#' Normalise data
#'
#' Normalises the read abundances to the total OTU counts (calculate percentages).
#'
#' @usage amp_normalise(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' 
#' @export
#' 
#' @author Rasmus Kirkegaard


amp_normalise <- function(data) {
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
  
  # Normalise data to percent
  data$abund <- apply(data$abund,2, function(x) 100*x/sum(x)) %>% as.data.frame() 

  return(data)
}
