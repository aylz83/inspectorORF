.observe_example_data_events <- function(input, output, session, files, tracks)
{
  ### RAW example data
  shiny::observeEvent(input$raw_example_data,
  {
    example_rna <- system.file("example_data", "control_rna_tracks.bed.gz", package = "inspectorORF")
    rna_name <- "rna_reads"

    example_ribo <- system.file("example_data", "control_psites_for_ORFquant", package = "inspectorORF")
    ribo_name <- "p_sites"

    files$rna_reads <- example_rna
    files$p_sites <- example_ribo

    new_rna_row <- data.frame(
      Name = rna_name,
      File = basename(example_rna),
      stringsAsFactors = FALSE
    )
    
    new_ribo_row <- data.frame(
      Name = ribo_name,
      File = basename(example_ribo),
      stringsAsFactors = FALSE
    )

    rnaseq_data$table <- new_rna_row
    riboseq_data$table <- new_ribo_row

    files$gtf_file <- system.file("example_data", "annotation_subset.gtf", package = "inspectorORF")
    files$two_bit <- system.file("example_data", "chr12.2bit", package = "inspectorORF")

    shiny::updateTextAreaInput(
      session,
      inputId = "id_input",
      value = "ENST00000343702.9"
    )
  })

  ### Processed example data
  shiny::observeEvent(input$bed_example_data,
  {
    example_data <- system.file("example_data", "example_tracks.bed", package = "inspectorORF")

    files$gtf_file <- system.file("example_data", "annotation_subset.gtf", package = "inspectorORF")
    files$two_bit <- system.file("example_data", "chr12.2bit", package = "inspectorORF")

    tracks$inspector_object <- inspectorORF::import_transcript_bed(
    	example_data,
    	gtf_file = files$gtf_file,
    	genome_file = files$two_bit
    )

    shinydashboard::updateTabItems(session, "tabs", selected = "data")
  })
}
