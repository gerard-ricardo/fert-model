launchMyApp <- function() {
  appDir <- system.file("custom_settings", package = "Fertmod")
  if (appDir == "") {
    stop("Shiny app not found in the package")
  }
  shiny::runApp(appDir)
}
