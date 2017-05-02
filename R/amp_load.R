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
  #Must be data frames
  otutable <- as.data.frame(otutable)
  metadata <- as.data.frame(metadata)
  colnames(metadata)[1] <- "SeqID"
  
  # Remove whitespace from the otutable as this will break the structure of the taxonomy
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  otutable$Kingdom<-trim(as.character(otutable$Kingdom))
  otutable$Phylum<-trim(as.character(otutable$Phylum))
  otutable$Class<-trim(as.character(otutable$Class))
  otutable$Order<-trim(as.character(otutable$Order))
  otutable$Family<-trim(as.character(otutable$Family))
  otutable$Genus<-trim(as.character(otutable$Genus))
  otutable$Species<-trim(as.character(otutable$Species))
  
  #metadata: order rows by rownames
  metadata = suppressWarnings(as.data.frame(as.matrix(metadata)))
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[order(rownames(metadata)), ]
  
  #Only alphanumeric characters in metadata column names, replace others with _
  colnames(metadata) <- str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")
  
  #abundance: all columns from otutable except the first and last 7 and convert to numeric for downstream compliance
  abund <- lapply(otutable[,2:(ncol(otutable) - 7)], as.numeric) %>% as.data.frame(check.names = F)
  rownames(abund) <- otutable$OTU # Add OTU as rownames

  #tax: the last 7 columns from otutable to factor, order rows by rownames and order columns by taxonomic rank(not alphabetically)
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] 
                    ,OTU = otutable$OTU, row.names = otutable$OTU)
  
  #data: return the data in a combined list w or w/o refseq. Load
  if(!is.null(refseq) & class(refseq) == "DNAStringSet") {
    data <- list(abund = abund, tax = tax, metadata = metadata, refseq = refseq)
  } else if(!is.null(refseq) & !class(refseq) == "DNAStringSet") {
    stop("The reference sequences must be loaded with readDNAStringSet() from the biostrings package.")
  } else if(is.null(refseq)) {
    data <- list(abund = abund, tax = tax, metadata = metadata)
  }
  
  #check if metadata and otutable match
  if(!all(rownames(data$metadata) %in% colnames(data$abund))) {
    stop("The sample names in metadata do not match those in otutable")
  }
  return(data)
}