.observe_process_events <- function(input, output, session, files, tracks)
{
  ### PROCESSING
  shiny::observeEvent(input$process_tracks,
  {
    tracks$filter_ids <- strsplit(input$id_input, "\n")[[1]]
    
    tracks$merged_tracks <- inspectorORF::merge_RNA_tracks_with_ORFquant(
    	rna_reads = files$rna_reads,
    	orfquant_psites = files$p_sites,
      extra = files$extra_reads,
      extra_orfquant_psites = files$extra_orfquant
    )

    tracks$inspector_object <- inspectorORF::get_id_tracks(
    	tracks$merged_tracks,
    	gtf_file = files$gtf_file,
    	genome_file = files$two_bit,
    	transcript_ids = tracks$filter_ids,
    	framed_tracks = tracks$framed_tracks
    )

    shinyjs::enable("export_tracks")
  })
}
