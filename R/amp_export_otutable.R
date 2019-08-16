#' Export OTU-table
#'
#' Export otutable (including taxonomy) from an ampvis2 object as CSV using \code{\link{write.table}}.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param filename File name of the exported OTU-table WITHOUT extension. (\emph{default:} \code{"exported_otutable"})
#' @param extension File extension without \code{"."}. (\emph{default:} \code{"csv"})
#' @param md5 (\emph{logical}) Compute md5 sum of the data (the whole object, not just otutable) and append to the filename. (\emph{default:} \code{FALSE})
#' @param sep Separator passed directly to \code{\link{write.table}}. (\emph{default:} \code{"\t"})
#' @param id Name the samples using a variable in the metadata.
#' @param sort_samples Vector to sort the samples by.
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{FALSE})
#' @param ... Additional arguments passed to \code{\link{write.table}} other than \code{sep} and \code{row.names}.
#'
#' @export
#'
#' @importFrom dplyr arrange mutate select desc everything
#' @importFrom digest digest
#'
#' @examples
#' # Load example data
#' data("AalborgWWTPs")
#'
#' # Export OTU-table
#' \dontrun{
#' amp_export_otutable(AalborgWWTPs, md5 = TRUE, filename = "AalborgWWTPs_otutable", sep = "\t")
#' }
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export_otutable <- function(data,
                                filename = "exported_otutable",
                                extension = "csv",
                                md5 = FALSE,
                                sep = "\t",
                                id = NULL,
                                sort_samples = NULL,
                                normalise = FALSE,
                                ...) {

  ### Data must be in ampvis2 format
  is_ampvis2(data)

  if (isTRUE(normalise)) {
    data <- normaliseTo100(data)
  }

  if (!is.null(id)) {
    ## Test if the ID exists in the metadata
    if (!(id %in% colnames(data$metadata))) {
      ametadata <- paste(colnames(data$metadata), collapse = ", ")
      stop(paste(id, "not found in metadata.\n\nAvailable metadata is: ", ametadata), call. = FALSE)
    }

    ## Test if the ID is unique for each sample
    if (length(unique(data$metadata[, id])) != length(colnames(data$abund))) {
      stop(paste(id, "is not unique for each sample"), call. = FALSE)
    }

    ## Re-arrange after coloumns after metadata
    re <- as.character(data$metadata[, 1])
    data$abund <- data$abund[, re]

    ## Add new sample names
    colnames(data$abund) <- as.character(unlist(data$metadata[, id]))
  }

  if (!is.null(sort_samples)) {

    ## Test if the ID is unique for each sample
    if (length(sort_samples) != length(colnames(data$abund))) {
      stop(paste("`sort_samples` does not match `id`"), call. = FALSE)
    }

    data$abund <- data$abund[, sort_samples]
  }

  # merge abundances and taxonomy by rownames
  e_bak <- merge(data$abund, data$tax, by = "row.names", all = TRUE, sort = FALSE)

  # remove first column (row.names) and order by OTU read counts across all samples
  e_bak2 <- e_bak %>%
    select(-1) %>%
    mutate(sum = rowSums(e_bak[, colnames(data$abund), drop = FALSE])) %>%
    arrange(desc(sum)) %>%
    select(-sum)

  # Append md5 sum to the filename just before the extenstion. Fx "../exported_otutable" will result in ../exported_otutable_md5sum.csv
  write.table(select(e_bak2, OTU, everything()),
    file = ifelse(md5,
      sprintf("%s_%s.%s", filename, digest::digest(data), extension),
      sprintf("%s.%s", filename, extension)
    ),
    quote = F,
    row.names = F,
    sep = sep,
    ...
  )
}
