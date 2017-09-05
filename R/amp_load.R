#' Load data for ampvis functions
#'
#' This function reads an OTU-table and corresponding sample metadata, and returns a list for use in all ampvis functions. It is therefore required to load data with \code{amp_load()} before any other ampvis functions can be used.
#'
#' @usage amp_load(otutable = dataframe, metadata = dataframe)
#'
#' @param otutable (\emph{required}) An OTU-table (data frame), where the last 7 rows is the taxonomy.
#' @param metadata (\emph{required}) A metadata (dataframe) with information about the samples.
#' @param refseq Reference sequences for all OTUs as loaded with \code{\link[Biostrings]{readDNAStringSet}} from the \code{\link{Biostrings}} bioconductor package.
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @import Biostrings
#' @import stringr
#' @export
#' 
#' @details The \code{amp_load()} function validates and corrects the provided data frames in different ways to make it suitable for the rest of the ampvis functions. It is important that the provided data frames match the requirements as described in the following sections to work properly.
#' 
#' @section The OTU-table:
#' The OTU-table contains information about the OTUs, their assigned taxonomy and their read counts in each sample. The provided OTU-table must be a data frame with the following requirements:
#' 
#' \itemize{
#'   \item The rows are OTU IDs and the columns are samples.
#'   \item The last 7 columns are the corresponding taxonomy assigned to the OTUs, named \code{"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}.
#'   \item The rownames of the data frame are the OTU IDs. Thus, the first column should be a sample, NOT the OTU IDs.
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
#' \code{amp_load()} will automatically use the sample IDs in the first column as rownames, but it is important to also have an actual column with sample IDs, so it is possible to fx group by that column during analysis.
#' 
#' A minimal example is available with \code{data("example_metadata")}.
#' 
#' @examples 
#' #Be sure to use the correct function to load your .csv files, see ?read.table()
#' 
#' \dontrun{
#' #Read the OTU-table as a data frame, top row will be used as column names (header = TRUE) 
#' #and first column as rownames. Taxonomy often has blank cells in the Species column, fill = TRUE
#' #makes sure it gets treated as such.
#' myotutable <- read.delim("data/otutable.txt", 
#'                         header = TRUE,
#'                         fill = TRUE,
#'                         stringsAsFactors = FALSE,
#'                         check.names = FALSE,
#'                         row.names = 1) %>% as.data.frame() 
#'                         
#' #Read the metadata, often an excel sheet. If .csv make sure the first column will be kept and NOT 
#' #loaded as rownames! The top row should be loaded column names
#' mymetadata <- read_excel("data/metadata.xlsx",
#'                          col_names = TRUE) %>% as.data.frame() 
#'                          
#' #Combine the data with amp_load() to make it compatible with ampvis2 functions.
#' #Uncomment the refseq line to load reference sequences (not required).
#' d <- amp_load(otutable = myotutable,
#'               metadata = mymetadata
#'               #, refseq = readDNAStringSet("data/otus.fa", format = "fasta") 
#'               )
#'               
#' #Show a short summary about the data by simply typing the name of the object in the console
#' d
#' }
#' 
#' #Minimal example otutable:
#' data("example_otutable")
#' example_otutable
#' 
#' #Minimal example metadata:
#' data("example_metadata")
#' example_metadata
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, refseq = NULL){
  ###check data
  otutable <- as.data.frame(otutable)
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata[,1]
  
  ### Check if refseq data is in the right format
  if(!is.null(refseq) & !class(refseq) == "DNAStringSet") {
    stop("The reference sequences must be loaded with readDNAStringSet() from the Biostrings bioconductor package.")
  }
  
  ### First column in OTU-table should be a sample, NOT the OTU ID's
  if (rownames(otutable) == otutable[,1]) {
    stop("The rownames of the OTU-table should not be identical to the first column.")
  }
  
  if (rownames(otutable) == c(1:nrow(otutable))) {
    stop("The rownames of the OTU-table do not seem to be OTU ID's:\n", as.character(head(rownames(otutable))))
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
  abund <- lapply(otutable[,1:(ncol(otutable) - 7)], as.numeric) %>% as.data.frame(check.names = FALSE, row.names = rownames(otutable))

  ### check if metadata and otutable match
  if(!all(rownames(metadata) %in% colnames(abund)) | nrow(metadata) != ncol(abund)) {
    stop("The samples in the metadata do not match those in the otutable.\nCheck that you have loaded matching files and that they meet the requirements described in ?amp_load().\nNumber of samples in metadata: ", nrow(metadata), "\nNumber of samples in otutable: ", ncol(otutable[,1:(ncol(otutable) - 7)]))
  }
  
  ### Abundance: re-arrange columns in the same order as the metadata
  abund <- abund[,as.character(metadata[,1])]
  
  ### tax: the last 7 columns from otutable
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] 
                    ,OTU = rownames(otutable)) #with added OTU column
  tax <- tax[order(rownames(tax)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  
  ###data: return the data in a combined list w or w/o refseq.
  if(!is.null(refseq)) {
    data <- list(abund = abund, tax = tax, metadata = metadata, refseq = refseq)
  } else {
    data <- list(abund = abund, tax = tax, metadata = metadata)
  }
  
  class(data) <- "ampvis2" #Our own "ampvis2" format yay!
  return(data)
}