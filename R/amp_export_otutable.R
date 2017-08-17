#' Export otutable from an ampvis object
#'
#' Export otutable from an ampcis object.
#'
#' @usage amp_export_otutable(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
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
  metadata <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  if(!is.null(id)){
    
    ## Test if the ID exists in the metadata
    if( !(id %in% colnames(metadata)) ){
      ametadata <- paste(colnames(metadata), collapse = ", ")
      stop(paste(id, "not found in metadata.\n\nAvailable metadata is: ", ametadata))
    } 
    
    ## Test if the ID is unique for each sample
    if( length(unique(metadata[,id])) != length(colnames(abund)) ){
      stop(paste(id, "is not unique for each sample"))
    } 

    ## Re-arrange after coloumns after metadata
    re <- as.character(metadata[,1])
    abund <- abund[,re]
    
    ## Add new sample names
    colnames(abund) <- as.character(unlist(data[["metadata"]][,id]))
  }
  
  if(!is.null(sort.samples)){
 
    ## Test if the ID is unique for each sample
    if( length(sort.samples) != length(colnames(abund)) ){
      stop(paste("`sort.samples` does not match `id`"))
    } 
      
    abund <- abund[,sort.samples]
  }
  
  e_bak <- cbind.data.frame(abund, tax)
  
  e_bak2 <- mutate(e_bak, 
                   sum = rowSums(e_bak[,1:nrow(data[["metadata"]])])) %>%
    arrange(desc(sum)) %>%
    select(-sum)
  
  write.table(e_bak2, file = file, quote = F, row.names = F, sep = "\t")
}
