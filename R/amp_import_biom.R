#' @title (Defunct) Import OTU table from BIOM file
#' @description Since ampvis2 v2.7.6 \code{amp_import_biom()} has been implemented directly in \code{\link{amp_load}} instead. Don't use.
#' @description Reads a .biom file and converts it into an OTU table compatible with ampvis2. This enables support for both the \href{http://qiime.org/}{QIIME} and \href{https://www.mothur.org/}{mothur} bioinformatic pipelines, as both software tools can output data in the \href{http://biom-format.org/}{BIOM format} (for mothur see \href{https://www.mothur.org/wiki/Make.biom}{make.biom}). Utilises the \href{https://github.com/joey711/biomformat}{biomformat} package, so both the JSON and HDF5 versions of the BIOM format are supported.
#'
#' @param ... Capture and ignore all args
#' 
#' @return A data frame
#' @export
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
amp_import_biom <- function(...) {
  args <- list(...)
  .Defunct(new = "amp_load(otutable = \"path/to/file.biom\")")
}
