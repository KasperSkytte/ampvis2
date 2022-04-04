#' Merge ampvis2 object(s)
#'
#' @description Merge any number of ampvis2 objects into a single object.
#'
#' @param ... (required) Any number of ampvis2-class objects to merge
#' @param by_refseq (recommended) Merge by exact matches between DNA reference sequences. The full DNA sequences will then be used as the new names in the output. (\emph{default:} \code{TRUE})
#'
#' @return An ampvis2-class object
#' @importFrom data.table rbindlist
#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#'
#' @details
#' It's important to ensure that the taxonomy for all OTU's across data sets is generated in the exact same way with the same database.
#' When \code{by_refseq = FALSE} it's likewise important to ensure that OTU ID's are not arbitrary between data sets and that they are corresponding to the same sequences across data sets (objects).
#' When \code{by_refseq = TRUE} the full DNA sequences will be used as the new OTU ID's. If the length of the names is a problem you can manually adjust the names in the output object in the \code{data$tax$OTU} column as well as setting the row names on both the \code{data$abund} and \code{data$tax} data frames. They must all be identical.
#' @export
#'
#' @examples
#' data("MiDAS")
#'
#' # summary of samples from 2010-2012
#' amp_subset_samples(
#'   MiDAS,
#'   Year %in% c("2010", "2011", "2012")
#' )
#'
#' # now merge individual objects and verify summary is the same
#' d_2010 <- amp_subset_samples(
#'   MiDAS,
#'   Year %in% "2010"
#' )
#' d_2011 <- amp_subset_samples(
#'   MiDAS,
#'   Year %in% "2011"
#' )
#' d_2012 <- amp_subset_samples(
#'   MiDAS,
#'   Year %in% "2012"
#' )
#'
#' amp_merge_ampvis2(
#'   d_2010,
#'   d_2011,
#'   d_2012
#' )
amp_merge_ampvis2 <- function(...,
                              by_refseq = TRUE) {
  obj_list <- list(...)

  # all objects must be ampvis2-class objects
  if (!all(sapply(obj_list, inherits, "ampvis2"))) {
    stop("One or more objects is not an ampvis2-class object", call. = FALSE)
  }

  # all objects must be either normalised or not, cant have both
  normalised <- obj_list %>%
    lapply(
      attr,
      which = "normalised",
      exact = TRUE
    ) %>%
    unlist()

  if (sum(normalised) > 0L & sum(normalised) != length(obj_list)) {
    stop("All objects must be either normalised or not, not mixed", call. = FALSE)
  }

  # no duplicate samples between objects are allowed (check abund)
  obj_list %>%
    lapply(
      function(x) {
        colnames(x[["abund"]])
      }
    ) %>%
    unlist() %>%
    duplicated() %>%
    any() %>%
    if (.) {
      stop("One or more samples occurs more than once between the objects (according to abundance table)", call. = FALSE)
    }

  # no duplicate samples between objects are allowed (check sample metadata)
  obj_list %>%
    lapply(
      function(x) {
        x[["metadata"]][[1]]
      }
    ) %>%
    unlist() %>%
    duplicated() %>%
    any() %>%
    if (.) {
      stop("One or more samples occurs more than once between the objects (according to sample metadata)", call. = FALSE)
    }

  if (isTRUE(by_refseq)) {
    # ensure all objects have refseqs loaded
    has_refseq <- sapply(
      obj_list,
      function(obj) {
        inherits(obj$refseq, c("DNAbin", "AAbin"))
      }
    )
    if (!all(has_refseq)) {
      stop("All objects must have DNA sequences loaded to be \
        able to merge by DNA sequence (recommended). Otherwise \
        merge by OTU name by setting by_refseq = FALSE if you're \
        sure it makes sense for your data.", call. = FALSE)
    }

    obj_list <- obj_list %>%
      lapply(
        function(obj) {
          obj$abund$OTU <- rownames(obj$abund)
          refseq_chr <- obj$refseq %>%
            as.character() %>%
            lapply(paste, collapse = "") %>%
            unlist(use.names = TRUE)
          obj$abund <- obj$abund[names(refseq_chr), ]
          obj$abund$OTU <- refseq_chr -> rownames(obj$abund)

          obj$tax <- obj$tax[names(refseq_chr), ]
          obj$tax$OTU <- refseq_chr -> rownames(obj$tax)
          obj
        }
      )
  }

  # merge abundance tables
  abund <- obj_list %>%
    lapply(
      function(obj) {
        obj$abund$OTU <- rownames(obj$abund)
        obj$abund
      }
    ) %>%
    reduce(
      full_join,
      by = "OTU"
    )

  # merge metadata
  # ensure the first sample ID column are named the same, just use that of the first obj
  idcolname <- colnames(obj_list[[1]][["metadata"]])[1]
  metadata <- obj_list %>%
    lapply(
      function(obj) {
        colnames(obj[["metadata"]])[1] <- idcolname
        obj[["metadata"]]
      }
    ) %>%
    rbindlist(
      fill = TRUE
    )

  # merge taxonomy
  taxonomy <- obj_list %>%
    lapply(
      `[[`,
      "tax"
    ) %>%
    rbindlist(
      fill = TRUE
    ) %>%
    unique()

  # check for contradictions in taxonomy. If an OTU occurs more than once
  # after running unique() in the code above there must be conflicts where
  # where taxonomy is different
  if (any(duplicated(taxonomy[["OTU"]]))) {
    stop("Conflicting taxonomy between one or more OTU's across all objects. Have they been classified in the exact same way across all objects?")
  }

  # merge refseq
  fasta <- sapply(obj_list, function(obj) {
    strsplit(obj$tax$OTU, "")
  })

  # load and return
  amp_load(
    otutable = abund,
    taxonomy = taxonomy,
    metadata = metadata,
    fasta = NULL,
    tree = NULL,
    pruneSingletons = FALSE
  )
}
