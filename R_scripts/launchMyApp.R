launchMyApp <- function() {
  appDir <- system.file("R_scripts", package = "Fertmod")
  if (appDir == "") {
    stop("Shiny app not found in the package")
  }
  shiny::runApp(appDir)
}
