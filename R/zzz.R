#' .onAttach
#'
#' @param lib lib
#' @param pkg pkg
#'
#' @importFrom utils packageVersion
#' @importFrom data.table setDTthreads
#' @keywords internal
.onAttach <- function(lib, pkg) {
  options(scipen = 6)
  setDTthreads(0)
  # Check for new github release version. (Not master branch version, release version!)
  if (!interactive()) {
    return()
  } else {
    installed_version <- as.character(utils::packageVersion(pkg))
    packageStartupMessage(
      paste0(
        "This is ",
        pkg,
        " version ",
        installed_version
      ),
      appendLF = FALSE
    )
    gitHubUser <- "kasperskytte"
    tryCatch(
      {
        DESCRIPTION <- readLines(
          paste0(
            "https://raw.githubusercontent.com/",
            gitHubUser,
            "/",
            pkg,
            "/master/DESCRIPTION"
          )
        )
        remote_version <- gsub("Version:\\s*", "", DESCRIPTION[grep("Version:", DESCRIPTION)])
        if (installed_version < remote_version) {
          packageStartupMessage(
            paste0(" (newer version ", remote_version, " available)"),
            appendLF = FALSE
          )
        } else if (installed_version >= remote_version) {
          packageStartupMessage(" (up to date)", appendLF = FALSE)
        }
      },
      error = function(e) {
        warning("Can't reach GitHub to check for new version just now. Trying again next time.", call. = FALSE)
      },
      warning = function(e) {
        warning("Can't reach GitHub to check for new version just now. Trying again next time.", call. = FALSE)
      }
    )
    packageStartupMessage(". Great documentation is available at the ampvis2 website: https://kasperskytte.github.io/ampvis2/", appendLF = TRUE)
  }
}
