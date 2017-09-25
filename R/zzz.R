.onAttach <- function(lib, pkg)  {
  packageStartupMessage("This is ampvis version ",
                        utils::packageDescription("ampvis2",
                                                  fields="Version"),
                        appendLF = TRUE)
}
