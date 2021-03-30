#' @title Import OTU table from BIOM file
#' @description Reads a .biom file and converts it into an OTU table compatible with ampvis2. This enables support for both the \href{http://qiime.org/}{QIIME} and \href{https://www.mothur.org/}{mothur} bioinformatic pipelines, as both software tools can output data in the \href{http://biom-format.org/}{BIOM format} (for mothur see \href{https://www.mothur.org/wiki/Make.biom}{make.biom}). Utilises the \href{https://github.com/joey711/biomformat}{biomformat} package, so both the JSON and HDF5 versions of the BIOM format are supported.
#'
#' @param file Path to the .biom file.
#'
#' @return A data frame
#' @export
#' @importFrom dplyr bind_rows
#'
#' @seealso
#' \code{\link{amp_load}}, \code{\link{amp_import_usearch}}
#'
#' @examples
#' \dontrun{
#' # First import the BIOM format OTU table:
#' biom_otutable <- amp_import_biom("path/to/file.biom")
#'
#' # Then use amp_load() with or without metadata as normal:
#' d <- amp_load(biom_otutable, metadata)
#' }
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
amp_import_biom <- function(file) {
  checkReqPkg(
    "biomformat",
    " Please install with:\n  install.packages(\"BiocManager\"); BiocManager::install(\"biomformat\")"
  )
  
  # Detect the file type and read the file
  if (tolower(tools::file_ext(file)) == "biom") {
    x <- biomformat::read_biom(file)

    # error if no taxonomy found in the file
    x$rows %>%
      lapply(function(row) {
        row$metadata$taxonomy
      }) %>%
      unlist(use.names = FALSE) %>%
      is.null() %>%
      all() %>%
      if (.) {
        stop("Cannot find the taxonomy of one or more OTU's in the provided .biom file", call. = FALSE)
      }

    # extract OTU read counts
    abund <- biomformat::biom_data(x) %>%
      as.matrix(check.names = FALSE) %>%
      as.data.frame(check.names = FALSE)

    # extract the taxonomy
    taxlist <- lapply(x$rows, function(x) {
      x$metadata$taxonomy
    })

    names(taxlist) <- lapply(x$rows, function(x) {
      x$id
    })

    tax <- as.data.frame(t(as.data.frame(taxlist, check.names = FALSE, stringsAsFactors = FALSE)))

    taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    if (ncol(tax) > 7L) {
      warning("Taxonomy had more than 7 levels for at least one OTU, discarding all excess levels except the first 7, and assuming they are Kingdom->Species. You may want to assure that this is correct.", call. = FALSE)
      tax <- tax[, 1:7, drop = FALSE]
    }

    # rename taxonomic levels
    colnames(tax) <- taxLevels[1:ncol(tax)]

    if (ncol(tax) < 7L) {
      warning("Taxonomy had less than 7 levels for all OTU's (Kingdom->Species), filling with NA from Species level and up.", call. = FALSE)
      taxSkeleton <- data.frame(matrix(ncol = 7, nrow = 0))
      colnames(taxSkeleton) <- taxLevels
      tax <- dplyr::bind_rows(taxSkeleton, tax)
    }

    # combine abundances and taxonomy and return
    otutable <- cbind(abund, tax) # no need for merge()
    return(otutable)
  } else if (tolower(tools::file_ext(file)) != "biom") {
    stop("The provided file is not in .biom format", call. = FALSE)
  }
}
