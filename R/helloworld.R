.onAttach <- function(lib, pkg)  {
  packageStartupMessage("This is optimus ",
                        utils::packageDescription("optimus",
                                                  field = "Version"),
                        ", roll out!",
                        appendLF = TRUE)
}
