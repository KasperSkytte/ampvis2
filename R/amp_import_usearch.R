#' @title Import OTU table from USEARCH pipelines
#' @description Reads an OTU count table and a corresponding sintax taxonomy table from a \href{http://www.drive5.com/usearch/}{USEARCH} analysis \href{http://www.drive5.com/usearch/manual/ex_miseq.html}{pipeline}, and then converts and combines them into an OTU table compatible with ampvis2.
#'
#' @param otutab Path to the tabulated OTU count table as created by \code{usearch10 -otutab ... -otutabout} commands.
#' @param sintax Path to the tabulated sintax taxonomy table as created by \code{usearch10 -sintax ... -tabbedout} commands.
#'
#' @return A data frame
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr left_join
#' @importFrom stringr str_extract
#' @examples
#' \dontrun{
#' # First import the usearch format OTU table:
#' usearch_otutable <- amp_import_usearch(otutab = "/path/to/otutab.txt", sintax = "/path/to/sintax.txt")
#'
#' # Then use amp_load() with or without metadata as normal:
#' d <- amp_load(usearch_otutable, metadata)
#' }
amp_import_usearch <- function(otutab, sintax) {
  # Read otutable (=read counts per sample table)
  counts <- data.table::fread(otutab,
    sep = "\t",
    fill = TRUE,
    header = TRUE,
    data.table = TRUE
  )
  colnames(counts)[1] <- "OTU"

  # Read taxonomy
  tax <- data.table::fread(sintax,
    sep = "\t",
    fill = TRUE,
    header = FALSE,
    data.table = TRUE
  )[, c(1, 4)]
  colnames(tax) <- c("OTU", "tax")

  # Separate each taxonomic level into individual columns.
  # This has to be done separately as taxonomic levels can be blank
  # in between two other levels.
  tax[, "Kingdom" := gsub("[dk]:", "k__", stringr::str_extract(tax, "[dk]:[^,]*"))]
  tax[, "Phylum" := gsub("p:", "p__", stringr::str_extract(tax, "p:[^,]*"))]
  tax[, "Class" := gsub("c:", "c__", stringr::str_extract(tax, "c:[^,]*"))]
  tax[, "Order" := gsub("o:", "o__", stringr::str_extract(tax, "o:[^,]*"))]
  tax[, "Family" := gsub("f:", "f__", stringr::str_extract(tax, "f:[^,]*"))]
  tax[, "Genus" := gsub("g:", "g__", stringr::str_extract(tax, "g:[^,]*"))]
  tax[, "Species" := gsub("s:", "s__", stringr::str_extract(tax, "s:[^,]*"))]
  tax <- tax[, -2]
  # the below would be more concise, but only works if all levels has a value,
  # fx d:test,p:test,o:test,f:test,g:test,s:test is missing "c:class" because
  # of low bootstrap value, and this would cause the other levels to be skewed and assigned to the wrong levels:
  # sintax[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := data.table::tstrsplit(tax, ",", fixed = TRUE)]

  # Join the counts and taxonomy, keep only OTU's in counts table
  otutable <- dplyr::left_join(counts, tax, by = "OTU")
  return(otutable)
}
