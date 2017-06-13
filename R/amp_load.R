#' Load data and convert to a list.
#'
#' Load data and convert to a list.
#'
#' @usage amp_load(otutable, metadata)
#'
#' @param otutable (required) A OTU table. The first row should be OTU name and the last 7 rows taxonomy.
#' @param metadata (required) A metadata file with sample names in first column.
#' @param refseq Reference sequences for all OTUs. Must be loaded with readDNAStringSet() from the biostrings package.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, refseq = NULL){
  
  otutable <- as.data.frame(otutable)
  colnames(otutable)[1] <- "OTU"
  
  metadata <- as.data.frame(metadata)
  colnames(metadata)[1] <- "SeqID"
  
  ### Only alphanumeric characters in metadata column names, replace others with _
  colnames(metadata) <- str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")
  
  ### Check that a "name" is in the first column of the otutable
  
  if(is.factor(otutable[,1]) == FALSE & is.character(otutable[,1]) == FALSE){
    stop("First column in the OTU table must be a OTU name (not integer)")
  }
  
  ### Abundance: all columns from otutable except the first and last 7 and convert to numeric for downstream compliance
  abund <- lapply(otutable[,2:(ncol(otutable) - 7)], as.numeric) %>% as.data.frame(check.names = F)
  rownames(abund) <- otutable$OTU # Add OTU as rownames
  
  
  ### Check if metadata and abundance data match
  if(!all((metadata$SeqID) %in% colnames(abund))) {
    stop("The sample names in metadata do not match those in otutable")
  }
  
  ### Check if taxonomic data has the correct names
  tax.names <- colnames(otutable[, (ncol(otutable) - 6):ncol(otutable)])
  expected.tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  if(!all(tax.names %in% expected.tax)) {
    stop(paste("The sample names in metadata do not match those in otutable.\nOtutable:", 
               paste(tax.names, collapse = ", "), 
               "\nExpected:", 
               paste(expected.tax, collapse =", ")))
  }
  
  ### Abundance: re-arrange columns in the same order as the metadata
  abund <- abund[,as.character(metadata[,1])]

  ### Check if OTUs are in the right format
  if(!is.null(refseq) & !class(refseq) == "DNAStringSet") {
    stop("The reference sequences must be loaded with readDNAStringSet() from the biostrings package.")
  }
  
  ### Remove whitespace from the otutable as this will break the structure of the taxonomy
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  otutable$Kingdom<-trim(as.character(otutable$Kingdom))
  otutable$Phylum<-trim(as.character(otutable$Phylum))
  otutable$Class<-trim(as.character(otutable$Class))
  otutable$Order<-trim(as.character(otutable$Order))
  otutable$Family<-trim(as.character(otutable$Family))
  otutable$Genus<-trim(as.character(otutable$Genus))
  otutable$Species<-trim(as.character(otutable$Species))
  
  ### tax: the last 7 columns from otutable
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] 
                    ,OTU = otutable$OTU, row.names = otutable$OTU)
  
  #data: return the data in a combined list w or w/o refseq.
  if(!is.null(refseq)) {
    data <- list(abund = abund, tax = tax, metadata = metadata, refseq = refseq)
  } else {
    data <- list(abund = abund, tax = tax, metadata = metadata)
  }
  
  return(data)
}