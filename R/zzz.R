.onAttach <- function(...){
  msg <- paste("Welcome! :)")
  packageStartupMessage(msg)
}
