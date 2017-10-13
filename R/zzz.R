.onAttach <- function(lib, pkg)  {
  #Check for new github release version. (Not master branch version, release version!)
  if (!interactive()) {
    return()
  } else {
    local_version <- utils::packageVersion("ampvis2")
    packageStartupMessage("This is ", pkg, " version ", local_version, ". Great documentation is available at the ampvis2 website: https://madsalbertsen.github.io/ampvis2/\n", appendLF = TRUE)
    if(requireNamespace("remotes", quietly = TRUE)) {
      tryCatch({
        github_ref <- remotes:::github_resolve_ref(
          remotes::github_release(), 
          remotes:::parse_git_repo("madsalbertsen/ampvis2@*release"))$ref
        github_version <- package_version(gsub("v", "", github_ref))
        if(local_version < github_version) {
          packageStartupMessage(
            "New release of ", pkg, " (", github_version, ") is available! Install the latest release with: \nremotes::install_github(\"madsalbertsen/ampvis2@*release\")\n\nRead the release notes at: https://github.com/MadsAlbertsen/ampvis2/releases/tag/", github_version)
        }
      }, error=function(e) {
        packageStartupMessage("Can't reach GitHub to check for new releases just now. Trying again next time. \n")
      })
    }
  }
}
