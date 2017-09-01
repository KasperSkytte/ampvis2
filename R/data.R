#' MiDAS example data
#'
#' The microbial community composition of the activated sludge of 55 different Danish Wastewater Treatment Plants from 2006-2013 as determined by 16S rRNA amplicon sequencing using the Illumina platform.
#'
#' @name MiDAS
#' @docType data
#'
#' @usage data(MiDAS)
#'
#' @format An object of class \code{"ampvis2"}. See \code{\link[ampvis2]{amp_load}}.
#'
#' @keywords data
#'
#' @references McIlroy et al. (2015); MiDAS: the field guide to the microbes of activated sludge, Database, Volume 2015, bav062
#' (\href{https://doi.org/10.1093/database/bav062}{https://doi.org/10.1093/database/bav062})
#'
#' @source \href{http://midasfieldguide.org}{MiDAS: Field Guide to the Microbes of Activated Sludge and Anaerobic Digesters}
#'
#' @examples
#' data(MiDAS)
#' MiDAS
NULL

#' A subset of the MiDAS example data
#'
#' A smaller subset of \code{\link{MiDAS}} with 50 samples from 2 Danish Wastewater Treatment Plants; Aalborg West and Aalborg East.
#'
#' @name AalborgWWTPs
#' @docType data
#'
#' @usage data(AalborgWWTPs)
#'
#' @format An object of class \code{"ampvis2"}. See \code{\link[ampvis2]{amp_load}}.
#'
#' @keywords data
#' 
#' @examples
#' data(AalborgWWTPs)
#' AalborgWWTPs
NULL

#' Functional information of microbes in the MiDAS database
#'
#' A data frame with known/unknown functional characteristics of the Genera in the MiDAS database. 
#'
#' @name MiF
#' @docType MiF
#'
#' @usage data(MiF)
#'
#' @format A data frame.
#' \describe{
#'   \item{POS}{Positive.}
#'   \item{NEG}{Negative.}
#'   \item{VAR}{Variable within species.}
#'   \item{NT}{Not tested.}
#' }
#'
#' @keywords data
#'
NULL