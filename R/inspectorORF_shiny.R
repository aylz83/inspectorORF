#' Launch an interactive application of inspectorORF
#'
#' Build and optionally launch the main application
#' @param launch Run immediately (default) or just return the shiny app object
#' @param host Host to start the shiny app on (defaults to 127.0.0.1)
#' @param port Port the shiny app will listen on (defaults to 4739)
#' @param launch_browser launch the browser automatically if TRUE (default)
#' @return a shiny object if launch = FALSE, otherwise, runs the app
#' @export
run_inspectorORF <- function(launch = TRUE, host = "127.0.0.1", port = 4739, launch_browser = TRUE)
{
  required_pkgs <- c("shiny", "shinyjs", "shinyFiles", "shinycssloaders", "shinydashboard", "reactable")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0)
  {
    stop(
      "The following required packages are missing: ", 
      paste(missing_pkgs, collapse = ", "), "\n",
      "Please install them with:\n",
      "install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))"
    )
  }

  app <- shiny::shinyApp(ui = .shiny_ui(), server = .shiny_server())

  if (launch)
  {
    shiny::runApp(
      app,
      launch.browser = launch_browser,
      host = host,
      port = port
    )
  }
  else
  {
    return(app)
  }
}


