.onAttach <- function(...){
  msg <- paste("Warning! By default model outputs are saved within the library",
               "directory. If this library is reinstalled, these files will be",
               "deleted.")
  packageStartupMessage(msg)
}
