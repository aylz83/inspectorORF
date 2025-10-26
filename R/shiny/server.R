.shiny_server <- function()
{
  function(input, output, session)
  {
    shinyjs::disable("export_tracks")
    options(shiny.maxRequestSize = 2048 * 1024^2)

    files <- shiny::reactiveValues(
      gtf = NULL,
      two_bit = NULL,
      rna_reads = NULL,
      p_sites = NULL,
      extra_reads = list(),
      extra_orfquant = list()
    )
  
    tracks <- shiny::reactiveValues(
      filter_ids = NULL,
      tx_object = NULL,
      framed_tracks = "p_sites"
    )

    .observe_upload_events(input, output, session, files, tracks)
    .observe_example_data_events(input, output, session, files, tracks)
    .observe_process_events(input, output, session, files, tracks)

    ### EXPORTING
    output$export_tracks <- shiny::downloadHandler(
      filename = function()
      {
        paste0("tracks_", Sys.Date(), ".bed")
      },
      content = function(file)
      {
        inspectorORF::save_tracks(tracks$inspector_object, file)
      }
    )
  }
}
