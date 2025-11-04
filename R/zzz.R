.onLoad <- function(libname, pkgname)
{
  loadNamespace("IRanges")
  loadNamespace("GenomicRanges")

  # Load nested R files manually
  shiny_dir <- system.file("R/shiny", package = pkgname)
  
  if (dir.exists(shiny_dir)) {
    r_files <- list.files(
      shiny_dir, 
      pattern = "\\.[rR]$", 
      recursive = TRUE, 
      full.names = TRUE
    )
    sapply(r_files, sys.source, envir = parent.env(environment()))
  }
}

