#' MiDAS example data
#'
#' The microbial community composition of the activated sludge of 55 different Danish Wastewater Treatment Plants from 2006-2013 as determined by 16S rRNA amplicon sequencing using the Illumina platform.
#'
#' @name MiDAS
#' @docType data
#'
#' @format An object of class \code{"ampvis2"}. See \code{\link[ampvis2]{amp_load}}.
#'
#' @keywords datasets
#'
#' @references McIlroy et al. (2015); MiDAS: the field guide to the microbes of activated sludge, Database, Volume 2015, bav062
#' (\doi{10.1093/database/bav062})
#'
#' @source \href{http://midasfieldguide.org}{MiDAS: Field Guide to the Microbes of Activated Sludge and Anaerobic Digesters}
#'
#' @examples
#' \dontrun{
#' data("MiDAS")
#' View(MiDAS)
#' }
NULL

#' A subset of the MiDAS example data
#'
#' A smaller subset of \code{\link{MiDAS}} with 50 samples from 2 Danish Wastewater Treatment Plants; Aalborg West and Aalborg East.
#'
#' @name AalborgWWTPs
#' @docType data
#'
#' @format An object of class \code{"ampvis2"}. See \code{\link[ampvis2]{amp_load}}.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data("AalborgWWTPs")
#' View(AalborgWWTPs)
#' }
NULL

#' A minimal example of an OTU-table
#'
#' A minimal example of an OTU-table that is compatible with \code{\link{amp_load}} and ampvis2 functions.
#'
#' @name example_otutable
#' @docType data
#'
#' @format a data frame.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data(example_otutable)
#' View(example_otutable)
#' }
NULL

#' A minimal example of sample metadata
#'
#' A minimal example of a metadata sheet that is compatible with \code{\link{amp_load}} and ampvis2 functions.
#'
#' @name example_metadata
#' @docType data
#'
#' @format a data frame.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data(example_metadata)
#' View(example_metadata)
#' }
NULL

#' A minimal example of taxonomy
#'
#' A minimal example of a taxonomy table that is compatible with \code{\link{amp_load}} and ampvis2 functions.
#'
#' @name example_taxonomy
#' @docType data
#'
#' @format a data frame.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data(example_taxonomy)
#' View(example_taxonomy)
#' }
NULL

#' Functional information about microbes at Genus level
#'
#' Grabbed from midasfieldguide.org at 2020-12-01
#'
#' @name midasfunctions_20201201
#' @docType data
#'
#' @format a data frame.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' data(midasfunctions_20201201)
#' View(midasfunctions_20201201)
#' }
NULL
