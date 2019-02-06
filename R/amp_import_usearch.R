#' @title Import OTU table from USEARCH pipelines
#' @description Reads an OTU count table and a corresponding sintax taxonomy table from a \href{http://www.drive5.com/usearch/}{USEARCH} analysis \href{http://www.drive5.com/usearch/manual/ex_miseq.html}{pipeline}, and then converts and combines them into an OTU table compatible with ampvis2.
#'
#' @param otutab Path to the tabulated OTU count table as created by \code{usearch10 -otutab ... -otutabout} commands.
#' @param sintax Path to the tabulated sintax taxonomy table as created by \code{usearch10 -sintax ... -tabbedout} commands.
#'
#' @return A data frame
#' @export
#' @importFrom plyr ldply
#' @importFrom dplyr mutate mutate_all select starts_with
#' @importFrom stringr str_c str_replace_all str_split
#' @examples
#' \dontrun{
#' # First import the usearch format OTU table:
#' usearch_otutable <- amp_import_usearch(otutab = "/path/to/otutab.txt", sintax = "/path/to/sintax.txt")
#'
#' # Then use amp_load() with or without metadata as normal:
#' d <- amp_load(usearch_otutable, metadata)
#' }
amp_import_usearch <- function(otutab, sintax) {
  otutab <- readLines(otutab)
  otutab[1] <- stringr::str_replace_all(otutab[1], "#", "")
  otutab <- stringr::str_split(otutab, "\t", simplify = TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate_all(as.character) %>%
    `colnames<-`(as.character(t(.[1, ]))) %>%
    .[-1, ]
  colnames(otutab)[1] <- "OTU"

  sintax <- readLines(sintax)
  sintax <- lapply(sintax, function(x) {
    if (str_detect(x, "\\t$|\\-$|\\+$")) paste0(x, "k__") else return(x)
  })
  tax <- stringr::str_replace_all(sintax, ".* |.*\t", "") %>%
    stringr::str_replace_all(c(
      "k:|d:" = "k__",
      "p:" = "p__",
      "c:" = "c__",
      "o:" = "o__",
      "f:" = "f__",
      "g:" = "g__",
      "s:" = "s__",
      "\"" = ""
    )) %>%
    stringr::str_split(",") %>%
    {
      `names<-`(., stringr::str_c(stringr::str_replace_all(sintax, " .*|\t.*", "")))
    } %>%
    lapply(function(x) {
      `names<-`(x, substr(x, 0, 3))
    }) %>%
    plyr::ldply(rbind)
  if (all(c("d__", "k__") %in% colnames(tax))) {
    tax <- dplyr::select(tax, -dplyr::starts_with("d__"))
  }
  if (!any(tolower(colnames(tax)) == "k__")) {
    tax <- dplyr::mutate(tax, Kingdom = NA)
  }
  if (!any(tolower(colnames(tax)) == "p__")) {
    tax <- dplyr::mutate(tax, Phylum = NA)
  }
  if (!any(tolower(colnames(tax)) == "c__")) {
    tax <- dplyr::mutate(tax, Class = NA)
  }
  if (!any(tolower(colnames(tax)) == "o__")) {
    tax <- dplyr::mutate(tax, Order = NA)
  }
  if (!any(tolower(colnames(tax)) == "f__")) {
    tax <- dplyr::mutate(tax, Family = NA)
  }
  if (!any(tolower(colnames(tax)) == "g__")) {
    tax <- dplyr::mutate(tax, Genus = NA)
  }
  if (!any(tolower(colnames(tax)) == "s__")) {
    tax <- dplyr::mutate(tax, Species = NA)
  }
  tax <- dplyr::mutate_all(tax, as.character)
  colnames(tax) <- c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otutable <- merge(otutab, tax, by = "OTU")
  return(otutable)
}
