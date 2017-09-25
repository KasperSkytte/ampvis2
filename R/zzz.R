.onAttach <- function(lib, pkg)  {
  if (!interactive()) return()
  local_version <- utils::packageVersion("ampvis2")
  packageStartupMessage("This is ", pkg, " version ", local_version, appendLF = TRUE)
  
  if(requireNamespace('devtools', quietly = TRUE)) {
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
