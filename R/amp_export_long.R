#' Transform ampvis2 objects into a long-format data frame
#'
#' @description Transforms any ampvis2 object into a long data frame (\code{data.table}) to facilitate custom data analysis. Only the elements OTU counts (\code{data$abund}), taxonomy (\code{data$tax}), and metadata (\code{data$metadata}) are used, not phylogenetic tree or DNA sequences.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param metadata_vars A character vector of sample metadata variables to include from \code{data$metadata}, or \code{NULL}. The first column (sample IDs) will always be included. Default includes all.
#' @param tax_levels A character vector of taxonomic levels to include from \code{data$tax}, or \code{NULL}. The OTU column will always be included. Default includes all.
#'
#' @return A \code{data.table} in long format.
#' @export
#'
#' @importFrom data.table setDT melt setcolorder
#' @examples
#' # load minimal example data
#' d <- amp_load(example_otutable, example_metadata)
#'
#' # transform d into a long-format data frame
#' d_long <- amp_export_long(d, metadata_vars = "Date", tax_levels = c("OTU", "Genus"))
#'
#' # print the data frame (data.table)
#' d_long
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
amp_export_long <- function(data,
                            metadata_vars = colnames(data$metadata),
                            tax_levels = colnames(data$tax)) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # find name of first (samples) column
  samplesCol <- colnames(data$metadata)[[1]]

  # always include OTU column from taxonomy and first column (sample IDs) from metadata
  data$tax <- data$tax[, unique(c("OTU", tax_levels)), drop = FALSE]

  # merge tax and abund by rownames
  otutable <- base::merge(
    x = data$tax,
    y = data$abund,
    by = 0,
    all = TRUE
  )
  setDT(otutable)

  # melt by OTU
  long_data <- melt(
    data = otutable,
    id.vars = colnames(data$tax),
    measure.vars = colnames(data$abund),
    value.name = "count",
    variable.name = samplesCol,
    variable.factor = FALSE
  )

  # merge with chosen sample metadata variables
  if (!is.null(metadata_vars)) {
    data$metadata <- data$metadata[, unique(c(samplesCol, metadata_vars)), drop = FALSE]
    setDT(data$metadata)
    long_data <- data$metadata[long_data, on = samplesCol]
  }

  # reorder and return
  setcolorder(long_data, c(samplesCol, "OTU", "count"))
  return(long_data)
}
