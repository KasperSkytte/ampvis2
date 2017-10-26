#' Load data for ampvis functions
#'
#' This function reads an OTU-table and corresponding sample metadata, and returns a list for use in all ampvis functions. It is therefore required to load data with \code{\link{amp_load}} before any other ampvis functions can be used.
#'
#' @usage amp_load(otutable = dataframe, metadata = dataframe)
#'
#' @param otutable (\emph{required}) OTU-table (data frame) with the read counts of all OTU's, where the last 7 columns is the taxonomy (Kingdom -> Species).
#' @param metadata (\emph{required}) Sample metadata (dataframe) with information about the samples. If none provided, dummy metadata will be created.
#' @param fasta (\emph{optional}) Path to a FASTA file with reference sequences for all OTU's in the OTU-table. (\emph{default:} \code{NULL})
#' @param tree (\emph{optional}) Phylogenetic tree of class \code{"phylo"} as loaded with \code{\link[ape]{read.tree}}. (\emph{default:} \code{NULL})
#' 
#' @return A list of class \code{"ampvis2"} with 3 to 5 elements.
#' @import ape
#' @import stringr
#' @import dplyr
#' @export
#' 
#' @details The \code{\link{amp_load}} function validates and corrects the provided data frames in different ways to make it suitable for the rest of the ampvis functions. It is important that the provided data frames match the requirements as described in the following sections to work properly.
#' 
#' @section The OTU-table:
#' The OTU-table contains information about the OTUs, their assigned taxonomy and their read counts in each sample. The provided OTU-table must be a data frame with the following requirements:
#' 
#' \itemize{
#'   \item The rows are OTU IDs and the columns are samples.
#'   \item The last 7 columns are the corresponding taxonomy assigned to the OTUs, named \code{"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}.
#'   \item The OTU ID's are expected to be in eiher the rownames of the data frame or in a column called "OTU". Otherwise the function will stop with a message.
#'   \item The column names of the data frame are the sample IDs, exactly matching those in the metadata, (and the last 7 columns named \code{Kingdom} -> \code{Species}, of course). Thus, the first row should contain the read counts of one of the OTUs in each sample, NOT the sample IDs and taxonomy.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#' 
#' A minimal example is available with \code{data("example_otutable")}.
#' 
#' @section The metadata:
#' The metadata contains additional information about the samples, for example where each sample was taken, date, pH, treatment etc, which is used to compare and group the samples during analysis. The amount of information in the metadata is unlimited, it can contain any number of columns (variables), however there are a few requirements: 
#' 
#' \itemize{
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the OTU-table.
#'   \item Column classes matter, categorical variables should be loaded either \code{as.character()} or \code{as.factor()}, and continuous variables \code{as.numeric()}. See below.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#' 
#' If for example a column is named "Year" and the entries are simply entered as numbers (2011, 2012, 2013 etc), then R will automatically consider these as numerical values (\code{as.numeric()}) and therefore the column as a continuous variable, while it is a categorical variable and should be loaded \code{as.factor()} or \code{as.character()} instead. This has consequences for the analysis as R treats them differently. Therefore either use the \code{colClasses = } argument when loading a csv file or \code{col_types = } when loading an excel file, or manually adjust the column classes afterwards with fx \code{metadata$Year <- as.character(metadata$Year)}.
#' 
#' The \code{\link{amp_load}} function will automatically use the sample IDs in the first column as rownames, but it is important to also have an actual column with sample IDs, so it is possible to fx group by that column during analysis. Any unmatched samples between the otutable and metadata will be removed. 
#' 
#' A minimal example is available with \code{data("example_metadata")}.
#' 
#' @section Reference sequences:
#' A fasta file with the raw sequences can be loaded as well, which will then be available in the refseq element of the ampvis2 object. These sequences will not be used in any ampvis2 function other than the two subset functions \code{\link{amp_subset_samples}} and \code{\link{amp_subset_taxa}}, so that they can be exported with \code{\link{amp_export_fasta}}. The fasta file is loaded with the \code{\link[ape]{read.FASTA}} function from the \code{\link{ape}} package. 
#' 
#' @examples 
#' #Be sure to use the correct function to load your .csv files, see ?read.table()
#' 
#' \dontrun{
#' #Read the OTU-table as a data frame. It is important to set check.names = FALSE.
#' myotutable <- read.delim("data/otutable.txt", check.names = FALSE) 
#'                         
#' #Read the metadata, often an excel sheet. If .csv make sure the first column will be kept and NOT 
#' #loaded as rownames! The top row should be loaded column names
#' mymetadata <- read_excel("data/metadata.xlsx", col_names = TRUE) 
#'                          
#' #Combine the data with amp_load() to make it compatible with ampvis2 functions.
#' #Uncomment the fasta line to load reference sequences (not required).
#' d <- amp_load(otutable = myotutable,
#'               metadata = mymetadata,
#'               fasta = "path/to/fastafile.fa" #optional
#'               )
#'               
#' #Show a short summary about the data by simply typing the name of the object in the console
#' d
#' }
#' 
#' #Minimal example metadata:
#' data("example_metadata")
#' example_metadata
#' 
#' #Minimal example otutable:
#' data("example_otutable")
#' example_otutable
#' 
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata = NULL, fasta = NULL, tree = NULL) {
  ### check data
  otutable <- as.data.frame(otutable)
  #create dummy metadata if none provided
  if(!is.null(metadata)) {
    metadata <- as.data.frame(metadata)
  } else if (is.null(metadata)) {
    metadata <- data.frame("SampleID" = colnames(otutable[,1:(ncol(otutable) - 7), drop = FALSE]), check.names = FALSE,
                           "DummyVariable" = "All samples")
    warning("No sample metadata provided, creating dummy metadata.")
  }
  rownames(metadata) <- as.character(metadata[,1])
  
  #check tree
  if(!is.null(tree) & !class(tree) == "phylo") {
    stop("The provided phylogenetic tree must be of class \"phylo\" as loaded with the ape::read.tree() function.")
  }
  
  ### OTU-table, check for OTU ID's
  if(any(tolower(colnames(otutable)) == "otu")) {
    otucolid <- which(tolower(colnames(otutable)) == "otu")
    rownames(otutable) <- as.character(otutable[,otucolid])
    otutable <- otutable[,-otucolid]
  } else if (all(rownames(otutable) %in% c(1:nrow(otutable))) & !any(tolower(colnames(otutable)) == "otu")) {
    stop("Cannot find OTU ID's. Make sure they are provided as rownames or in a column called \"OTU\".")
  }
  
  ### Only alphanumeric characters in metadata column names, replace others with "_", they may cause problems with ggplot2 groupings etc
  colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")
  
  ### Check if taxonomic data has the correct names
  tax.names <- colnames(otutable[, (ncol(otutable) - 6):ncol(otutable)])
  expected.tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if(!all(tax.names %in% expected.tax)) {
    stop(paste("The last 7 columns in the OTU-table must be the taxonomy (Kingdom->Species) and named accordingly\nCurrent:", 
               paste(tax.names, collapse = ", "), 
               "\nExpected:", 
               paste(expected.tax, collapse =", ")))
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
  
  ### Abundance: all columns from otutable except the first and last 7 and convert to numeric for downstream compliance
  abund <- lapply(otutable[,1:(ncol(otutable) - 7), drop = FALSE], as.numeric) %>%
    as.data.frame(check.names = FALSE, row.names = rownames(otutable))
  abund0 <- abund
  metadata0 <- metadata
  
  ### check if metadata and otutable match
  if(!all(rownames(metadata) %in% colnames(abund)) | !all(colnames(abund) %in% rownames(metadata))) {
    if (!any(colnames(abund) %in% rownames(metadata))) {
      # No samples match
      stop("No sample names match between metadata and otutable. Check that you have loaded matching files and that they meet the requirements described in ?amp_load(). Remember to use check.names = FALSE when loading the files.")
    } else {
      #Find matching samples between otutable and metadata
      sharedSamples <- dplyr::intersect(colnames(abund), rownames(metadata))
      
      #subset both based on shared samples
      abund0 <- abund[,match(sharedSamples, colnames(abund)), drop = FALSE]
      abund0 <- abund0[rowSums(abund0) > 0,] #after subset, rows with all 0's may be introduced, remove
      metadata0 <- metadata[match(sharedSamples, rownames(metadata)),, drop = FALSE]
      
      #Vectors with unique sample names in either
      metadataUniques <- metadata[-which(rownames(metadata) %in% rownames(metadata0)),, drop = FALSE] %>% rownames()
      abundUniques <- abund[,-which(colnames(abund) %in% colnames(abund0)), drop = FALSE] %>% colnames()
      
      #print the removed samples
      warning("Only ", ncol(abund0), " of ", length(unique(c(rownames(metadata), colnames(abund)))), 
              " unique sample names match between metadata and otutable. The following unmatched samples have been removed:", 
              ifelse(length(metadataUniques) > 0, paste0("\nmetadata (", length(metadataUniques), "): \n\t\"", paste(metadataUniques, collapse = "\", \""), "\""), ""),
              ifelse(length(abundUniques) > 0, paste0("\notutable (", length(abundUniques), "): \n\t\"", paste(abundUniques, collapse = "\", \""), "\""), ""))
    }
  } 
  
  ### tax: the last 7 columns from otutable
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] 
                    ,OTU = rownames(otutable)) #with added OTU column
  tax <- tax[which(rownames(abund) %in% rownames(abund0)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  tax[is.na(tax)] <- ""  
  
  ### data: return the data in a combined list w or w/o refseq. The rows of tax are ordered by the rownames of abund, and the columns of abund are ordered by the metadata rownames
  data <- list(abund = abund0[,rownames(metadata0), drop = FALSE], tax = tax[rownames(abund0),, drop = FALSE], metadata = metadata0)
  
  #apend refseq if provided
  if(!is.null(fasta)) {
    data[["refseq"]] <- ape::read.FASTA(file = fasta)
  } 
  
  #apend phylogenetic tree if provided
  if(!is.null(tree)) {
    data[["tree"]] <- tree
  }
  
  class(data) <- "ampvis2" #Our own "ampvis2" class, yay!
  return(data)
}
