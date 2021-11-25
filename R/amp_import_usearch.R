#' @title (Defunct) Import OTU table from USEARCH pipelines
#' @description Since ampvis2 v2.7.6 \code{amp_import_usearch()} has been implemented directly in \code{\link{amp_load}} instead. Don't use.
#' @description Reads an OTU count table and a corresponding sintax taxonomy table from a \href{http://www.drive5.com/usearch/}{USEARCH} analysis \href{http://www.drive5.com/usearch/manual/ex_miseq.html}{pipeline}, and then converts and combines them into an OTU table compatible with ampvis2.
#'
#' @param ... Capture and ignore all args
#'
#' @return A data frame
#' @export
#' @examples
#' \dontrun{
#' # First import the usearch format OTU table:
#' usearch_otutable <- amp_import_usearch(
#'   otutab = "/path/to/otutab.txt",
#'   sintax = "/path/to/sintax.txt"
#' )
#'
#' # Then use amp_load() with or without metadata as normal:
#' d <- amp_load(usearch_otutable, metadata)
#' }
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
amp_import_usearch <- function(...) {
  args <- list(...)
  .Defunct(new = "amp_load(
  otutable = \"path/to/zotutable.tsv\", \
  taxonomy = \"path/to/taxonomy.sintax\"
)")
}
