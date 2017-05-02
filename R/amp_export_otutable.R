#' Export otutable from an ampvis object
#'
#' Export otutable from an ampcis object.
#'
#' @usage amp_export_otutable(data)
#'
#' @param data (required) A ampvis object.
#' @param file Name of the file containing the exported otutable.
#' @param id Name the samples using any variable in the metadata.
#' @param sort.samples Vector to sort the samples by.
#' @param raw Export raw input instead of converting to percentages (default: F) 
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export_otutable <- function(data, file = "exported_otutable.txt", id = NULL, sort.samples = NULL, raw = F){
  
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  if(!is.null(id)){
    colnames(abund) <- as.character(unlist(data[["metadata"]][,id]))
  }
  
  if(!is.null(sort.samples)){
    abund <- abund[,sort.samples]
  }
  
  e_bak <- cbind.data.frame(abund, tax)
  
  e_bak2 <- mutate(e_bak, 
                   sum = rowSums(e_bak[,1:nrow(data[["metadata"]])])) %>%
    arrange(desc(sum)) %>%
    select(-sum)
  
  write.table(e_bak2, file = file, quote = F, row.names = F, sep = "\t")
}
