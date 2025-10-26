.observe_upload_events <- function(input, output, session, files, tracks)
{
  # Observe fileInput uploads
  shiny::observeEvent(input$gtf_file,
  {
    files$gtf <- input$gtf_file$datapath
  })

  shiny::observeEvent(input$sequence_file,
  {
    files$two_bit <- input$sequence_file$datapath
  })

  ##### RNA SEQ imports
  rnaseq_data <- shiny::reactiveValues(table = data.frame(
    Name = character(),
    File = character(),
    stringsAsFactors = FALSE
  ))

  output$rnaseq_name <- shiny::renderUI(
  {
    if (nrow(rnaseq_data$table) == 0)
    {
      # Only 1 file => no input, fixed name
      shiny::tags$p("Sample name will be 'rna_reads'")
    }
    else
    {
      shiny::textInput("rnaseq_name", "Sample name", value = "")
    }
  })

  shiny::observeEvent(input$add_rnaseq,
  {
    shiny::req(input$rnaseq_file)  # Make sure a file is uploaded
  
    # Get current table
    df <- rnaseq_data$table

    sample_name <- if (nrow(df) == 0)
    {
      file$rna_reads <- input$rnaseq_file$datapath
      "rna_reads"
    }
    else
    {
      file$extra_reads[[input$rnaseq_name]] <- input$rnaseq_file$datapath
      input$rnaseq_name
    }

    # Append new row
    new_row <- data.frame(
      Name = sample_name,
      File = input$rnaseq_file$name,
      stringsAsFactors = FALSE
    )
  
    rnaseq_data$table <- rbind(df, new_row)
  })

  output$rnaseq_table <- DT::renderDT(
    {
      rnaseq_data$table
    },
    options = list(pageLength = 6, scrollX = TRUE, searching = FALSE, lengthChange = FALSE)
  )

  #### Ribo-Seq imports
  riboseq_data <- shiny::reactiveValues(table = data.frame(
    Name = character(),
    File = character(),
    stringsAsFactors = FALSE
  ))

  output$riboseq_name <- shiny::renderUI(
  {
    if (nrow(riboseq_data$table) == 0)
    {
      # Only 1 file => no input, fixed name
      shiny::tags$p("Sample name will be 'p_sites'")
    }
    else
    {
      shiny::textInput("riboseq_name", "Sample name", value = "")
    }
  })

  shiny::observeEvent(input$add_orfquant,
  {
    shiny::req(input$riboseq_file)  # Make sure a file is uploaded
  
    # Get current table
    df <- riboseq_data$table

    sample_name <- if (nrow(df) == 0)
    {
      file$p_sites <- input$riboseq_file$datapath
      "p_sites"
    }
    else
    {
      file$extra_orfquant[[input$riboseq_name]] <- input$riboseq_file$datapath
      tracks$framed_tracks <- c(tracks$framed_tracks, input$riboseq_name)
      input$riboseq_name
    }

    # Append new row
    new_row <- data.frame(
      Name = sample_name,
      File = input$riboseq_file$name,
      stringsAsFactors = FALSE
    )
  
    riboseq_data$table <- rbind(df, new_row)
  })

  output$orfquant_table <- DT::renderDT(
    {
      riboseq_data$table
    },
    options = list(pageLength = 6, scrollX = TRUE, searching = FALSE, lengthChange = FALSE)
  )
}
