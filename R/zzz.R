.onAttach <- function(lib, pkg)  {
  require("ggplot2", quietly = TRUE, warn.conflicts = FALSE) #load ggplot2 without warnings
  if (!interactive()) {
    return()
  } else {
    packageStartupMessage("This is ", pkg, " version ", utils::packageVersion("ampvis2"), appendLF = TRUE)
    if(requireNamespace('remotes', quietly = TRUE)) {
      tryCatch({
        github_ref <- remotes:::github_resolve_ref(
          remotes::github_release(), 
          remotes:::parse_git_repo("madsalbertsen/ampvis2@*release"))$ref
        github_version <- package_version(gsub("v", "", github_ref))
        if(local_version != github_version) {
          packageStartupMessage(
            'New version of ', pkg, ' (', github_version, ') is available! Install the latest release with \n\"devtools::install_github("madsalbertsen/ampvis2@*release").\"')
        }
      }, error=function(e) {
        packageStartupMessage("Can't reach GitHub to check for new releases just now. Trying again next time. ")
      })
    }
  }
}
