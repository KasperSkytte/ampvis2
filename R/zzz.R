.onAttach <- function(lib, pkg)  {
  suppressWarnings(suppressMessages(requireNamespace("ggplot2", quietly = TRUE))) #load ggplot2 without warnings
  if (!interactive()) {
    return()
  } else {
    local_version <- utils::packageVersion("ampvis2")
    packageStartupMessage("This is ", pkg, " version ", local_version, appendLF = TRUE)
    if(requireNamespace("remotes", quietly = TRUE)) {
      tryCatch({
        github_ref <- remotes:::github_resolve_ref(
          remotes::github_release(), 
          remotes:::parse_git_repo("madsalbertsen/ampvis2@*release"))$ref
        github_version <- package_version(gsub("v", "", github_ref))
        if(local_version != github_version) {
          packageStartupMessage(
            "New release of ", pkg, " (", github_version, ") is available! Install the latest release of ", pkg, " with \nremotes::install_github(\"madsalbertsen/ampvis2@*release\").")
        }
      }, error=function(e) {
        packageStartupMessage("Can't reach GitHub to check for new releases just now. Trying again next time. ")
      })
    }
  }
  invisible()
}
