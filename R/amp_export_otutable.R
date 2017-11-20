#' Export OTU-table
#'
#' Export otutable (including taxonomy) from an ampvis2 object as CSV using \code{\link{write.table}}.
#'
#' @usage amp_export_otutable(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param filename File name of the exported OTU-table WITHOUT extension. (\emph{default:} \code{"exported_otutable"})
#' @param md5 (\emph{logical}) Compute md5 sum of the data (the whole object, not just otutable) and append to the filename. (\emph{default:} \code{FALSE})
#' @param sep Separator passed directly to \code{\link{write.table}}. (\emph{default:} \code{"\t"})
#' @param id Name the samples using a variable in the metadata.
#' @param sort_samples Vector to sort the samples by.
#' @param raw (\emph{logical}) Use raw counts instead of percentages. (\emph{default:} \code{TRUE})
#' @param ... Additional arguments passed to \code{\link{write.table}} other than \code{sep} and \code{row.names}.
#' 
#' @export
#' @import dplyr
#' @import digest
#' @import stringr
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Export OTU-table
#' \dontrun{
#' amp_export_otutable(AalborgWWTPs, md5 = TRUE, filename = "AalborgWWTPs_otutable", sep = "\t")
#' }
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export_otutable <- function(data,
                                filename = "exported_otutable",
                                md5 = FALSE,
                                sep = "\t",
                                id = NULL, 
                                sort_samples = NULL, 
                                raw = TRUE,
                                ...){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
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
  
  if(!is.null(sort_samples)){
 
    ## Test if the ID is unique for each sample
    if( length(sort_samples) != length(colnames(abund)) ){
      stop(paste("`sort_samples` does not match `id`"))
    } 
      
    abund <- abund[,sort_samples]
  }
  
  e_bak <- cbind.data.frame(abund, tax)
  
  e_bak2 <- mutate(e_bak, 
                   sum = rowSums(e_bak[,1:nrow(data[["metadata"]])])) %>%
    arrange(desc(sum)) %>%
    select(-sum)
  
  #Append md5 sum to the filename just before the extenstion. Fx "../exported_otutable" will result in ../exported_otutable_md5sum.csv
  write.table(select(e_bak2, OTU, everything()), file = ifelse(md5, sprintf("%s_%s.csv", filename, digest::digest(data)), paste0(filename, ".csv")), quote = F, row.names = F, sep = sep, ...)
}
