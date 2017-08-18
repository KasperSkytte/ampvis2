#' Load data for ampvis functions
#'
#' This function reads an OTU-table and corresponding sample metadata, and returns a list for use in all ampvis functions. It is therefore required to load data with \code{amp_load()} before any other ampvis functions can be used.
#'
#' @usage amp_load(otutable = dataframe, metadata = dataframe)
#'
#' @param otutable (\emph{required}) An OTU-table (data frame), where the last 7 rows is the taxonomy.
#' @param metadata (\emph{required}) A metadata (dataframe) with information about the samples.
#' @param refseq Reference sequences for all OTUs as loaded with \code{\link[biostrings]{readDNAStringSet()}} from the \code{biostrings} bioconductor package.
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' 
#' @export
#' 
#' @details The \code{amp_load()} function validates and corrects the provided data in different ways to make it suitable for the rest of the ampvis functions. It is important that the provided data frames match the requirements to work properly as described in the following sections. 
#' 
#' @section The OTU-table:
#' The OTU-table contains information about the OTUs, their assigned taxonomy and their read counts in each sample. The provided OTU-table must be a data frame with the following requirements:
#' 
#' \itemize{
#'   \item The rows are OTU IDs and the columns are samples.
#'   \item The last 7 columns are the corresponding taxonomy assigned to the OTUs, named \code{"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}.
#'   \item The rownames of the data frame are the OTU IDs. Thus, the first column should be a sample, NOT the OTU IDs.
#'   \item The column names of the data frame are the sample IDs, exactly matching those in the metadata, (the last 7 columns named \code{Kingdom} -> \code{Species}, of course). Thus, the first row should contain the read counts of an OTU in each sample, NOT the sample IDs and taxonomy.
#'   \item Generally avoid special characters in row- and column names.
#' }
#' 
#' @section The metadata:
#' The metadata contains additional information about the samples, for example where the sample was taken, date, pH, treatment etc, which is used to compare and group the samples during analysis. The amount of information in the metadata is unlimited, it can contain any number of columns (variables), however there are a few requirements: 
#' 
#' \itemize{
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the OTU-table.
#'   \item Column classes matter, categorical variables should be loaded either \code{as.character()} or \code{as.factor()}, and continuous variables \code{as.numeric()}.
#'   \item Generally avoid special characters in row- and column names.
#' }
#' 
#' If for example a column is named "Year" and the entries are simply entered as numbers (2011, 2012, 2013 etc), then R will consider these as numerical values (\code{as.numeric()}) and therefore the column as a continuous variable, while it is a categorical variable and should be loaded \code{as.factor()} or \code{as.character()} instead. This has consequences for the analysis as R treats them differently. Therefore either use the \code{colClasses = } argument when loading a csv file or \code{col_types = } when loading an excel file, or manually adjust the column classes afterwards with fx \code{metadata$Year <- as.character(metadata$Year)}.
#' 
#' \code{amp_load()} will automatically use the first column as rownames, but it is important to also have an actual column with sample IDs, so it is possible to fx group by that column during analysis.
#' 
#' @examples 
#' myotutable <- read.csv2("data/otutable.csv", 
#'                         header = TRUE,
#'                         stringsAsFactors = FALSE,
#'                         check.names = FALSE,
#'                         row.names = 1) %>% as.data.frame() 
#' mymetadata <- read_excel("data/metadata.xlsx",
#'                          col_names = TRUE) %>% as.data.frame() 
#' d <- amp_load(otutable = myotutable,
#'               metadata = mymetadata
#'               #, refseq = readDNAStringSet("data/otus.fa", format = "fasta") #optional 
#'               )
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, refseq = NULL){
  #check data
  otutable <- as.data.frame(otutable)
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata[,1]
  
  ### Only alphanumeric characters in metadata column names, replace others with _
  colnames(metadata) <- str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")
  
  ### Check if refseq data is in the right format
  if(!is.null(refseq) & !class(refseq) == "DNAStringSet") {
    stop("The reference sequences must be loaded with readDNAStringSet() from the biostrings bioconductor package.")
  }
  
  ### Check if taxonomic data has the correct names
  tax.names <- colnames(otutable[, (ncol(otutable) - 6):ncol(otutable)])
  expected.tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if(!all(tax.names %in% expected.tax)) {
    stop(paste("The last 7 columns in the OTU-table must be the taxonomy (Kingdom->Species) and named accordingly\nCurrent:", 
               paste(tax.names, collapse = ", "), 
               "\nExpected:", 
               paste(expected.tax, collapse =", ")))
  }
  
  ### Abundance: all columns from otutable except the first and last 7 and convert to numeric for downstream compliance
  abund <- lapply(otutable[,1:(ncol(otutable) - 7)], as.numeric) %>% as.data.frame(check.names = FALSE, row.names = rownames(otutable))

  ### Abundance: re-arrange columns in the same order as the metadata
  abund <- abund[,as.character(metadata[,1])]
  
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
                    ,OTU = rownames(otutable)) #with added OTU column
  tax <- tax[order(rownames(tax)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  
  #data: return the data in a combined list w or w/o refseq.
  if(!is.null(refseq)) {
    data <- list(abund = abund, tax = tax, metadata = metadata, refseq = refseq)
  } else {
    data <- list(abund = abund, tax = tax, metadata = metadata)
  }
  
  #calculate basic statistics and useful information about the data, print it
  c("Samples" = as.character(ncol(data$abund)),
    "OTUs" = as.character(nrow(data$abund)),
    "Total#Reads" = as.character(sum(data$abund)),
    "Min#Reads" = as.character(min(colSums(data$abund))),
    "Avg#Reads" = as.character(round(mean(colSums(data$abund)), digits = 2)),
    "Max#Reads" = as.character(max(colSums(data$abund)))) %>% print.table()
  cat("\nAssigned taxonomy:\n")
  c("Kingdom" = paste(round(sum(nchar(data$tax$Kingdom) > 3) / nrow(data$abund), digits = 2) * 100, "%", sep = ""),
    "Phylum" = paste(round(sum(nchar(data$tax$Phylum) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = ""),
    "Class" = paste(round(sum(nchar(data$tax$Class) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = ""),
    "Order" = paste(round(sum(nchar(data$tax$Order) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = ""),
    "Family" = paste(round(sum(nchar(data$tax$Family) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = ""),
    "Genus" = paste(round(sum(nchar(data$tax$Genus) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = ""),
    "Species" = paste(round(sum(nchar(data$tax$Species) > 3) / nrow(data$abund) * 100, digits = 2), "%", sep = "")) %>% print.table()
  
  return(data)
}