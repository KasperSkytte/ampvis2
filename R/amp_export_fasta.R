#' Export raw DNA sequences
#'
#' Export sequences from an ampvis2 object with added taxonomy.
#'
#' @usage amp_export_fasta(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param filename File name of the exported FASTA file. (\emph{default:} \code{"exported_sequences.fa"})
#' @param tax (\emph{logical}) Add taxonomic strings to the output or not. (\emph{default:} \code{TRUE})
#' 
#' @import ape
#' @import stringr
#' @export
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Export fasta sequences as "exported_sequences.fa"
#' \dontrun{
#' amp_export_fasta(AalborgWWTPs)
#' }
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export_fasta <- function(data, 
                             filename = "exported_sequences.fa", 
                             tax = TRUE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ### Reference sequences must be there!
  if(is.null(data$refseq)) {
    stop("No \"refseq\" element in the provided data.")
  }
  
  ### Check if refseq data is in the right format
  if(!is.null(data$refseq) & !class(data$refseq) == "DNAbin") {
    stop("The refseq element is not of class \"DNAbin\". The reference sequences must be loaded with ape::read.dna().")
  }
  
  t <- data[["refseq"]]
  
  if (tax == TRUE){
    #dont want duplicate taxonomy, if the refseq loaded already contains taxonomy.
    if(any(stringr::str_detect(names(t), ";"))) {
      names(t) <- stringr::str_split(names(t), ";", simplify = TRUE)[,1]
    }
    tax <- as.data.frame(data[["tax"]][,1:7])
    tax <- tax[names(t),]
    tax_s <- data.frame(lapply(tax, as.character), stringsAsFactors=FALSE)
    df_args <- c(tax_s, sep="; ")
    tax_sf <- do.call(paste, df_args)
    
    names(t) <- paste(names(t), tax_sf, sep = "; ")
  }
  
  ape::write.dna(t, file = filename, format = "fasta")
}
