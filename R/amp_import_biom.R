#' @title Import OTU table from BIOM file
#' @description Reads a .biom file and converts it into an OTU table compatible with ampvis2. This enables support for both the \href{http://qiime.org/}{QIIME} and \href{https://www.mothur.org/}{mothur} bioinformatic pipelines, as both software tools can output data in the \href{http://biom-format.org/}{BIOM format} (for mothur see \href{https://www.mothur.org/wiki/Make.biom}{make.biom}). Utilises the \href{https://github.com/joey711/biomformat}{biomformat} package, so both the JSON and HDF5 versions of the BIOM format are supported. 
#'
#' @param file Path to the .biom file.
#'
#' @return A data frame
#' @export
#' @importFrom biomformat read_biom biom_data
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' #First import the BIOM format OTU table:
#' biom_otutable <- amp_import_biom("path/to/file.biom")
#' 
#' #Then use amp_load() with or without metadata as normal:
#' d <- amp_load(biom_otutable, metadata)
#' }
amp_import_biom <- function(file) {
  #Detect the file type and read the file
  if(tolower(tools::file_ext(file)) == "biom") {
    x <- biomformat::read_biom(file)
    
    #error if no taxonomy found in the file
    if (all(sapply(sapply(x$rows, function(x) {
      x$metadata
    }), is.null)))
      stop("Cannot find taxonomy in the provided .biom file", call. = FALSE)
    
    #extract OTU read counts
    abund <- biomformat::biom_data(x) %>%
      as.matrix(check.names = FALSE) %>%
      as.data.frame(check.names = FALSE)
    
    #extract the taxonomy
    taxlist <- lapply(x$rows, function(x) {
      x$metadata$taxonomy
    })
    names(taxlist) = lapply(x$rows, function(x) {
      x$id
    })
    tax <- t(as.data.frame(taxlist, check.names = FALSE, stringsAsFactors = FALSE))
    colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    otutable <- cbind(abund, tax)
    return(otutable)
  } else if(tolower(tools::file_ext(file)) != "biom") {
    stop("The provided file is not in .biom format", call. = FALSE)
  }
}
