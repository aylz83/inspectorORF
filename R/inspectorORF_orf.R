#' @slot orfs A list of ORFs and associated data such as transcript_id, start, stop, framing etc.
setClass("inspectorORF_orflist",
         slots = c(orfs = "list"),
         prototype = list(orfs = list()))

#' Coerce the ORF object into a data.frame
#'
#' Coerce the ORF object into a data.frame
#' @param x The object created from inspectorORF::find_orfs()
#' @param row.names unused
#' @param optional unused
#' @param ... Additional arguments passed to or from other methods.
#' @return A data.frame consisting of the ORFs
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orfs <- inspectorORF::find_orfs(
#'   tracks,
#'   "ENST00000343702.9"
#' )
#' df <- test_orfs |> as.data.frame()
#' @export
as.data.frame.inspectorORF_orflist <- function(x, row.names = NULL, optional = FALSE, ...)
{
  orfs_df <- do.call(rbind, lapply(x@orfs, as.data.frame))

  # Check for unnamed columns
  if (is.null(colnames(orfs_df)) || any(colnames(orfs_df) == ""))
  {
    stop("Some or all columns are unnamed. Cannot convert to data.frame safely.")
  }

  orfs_df <- as.data.frame(orfs_df)
  orfs_df$aa_length <- orfs_df$nt_length / 3

  orfs_df
}

#' Print obtained ORFs
#'
#' Prints out the obtained ORFs as a data.frame
#' @param x The object created from inspectorORF::find_orfs()
#' @param ... Additional arguments passed to or from other methods.
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orfs <- inspectorORF::find_orfs(
#'   tracks,
#'   "ENST00000343702.9"
#' )
#' print(test_orfs)
#' @export
print.inspectorORF_orflist <- function(x, ...)
{
  print(
    as.data.frame(x),
    ...
  )
}

#' Get a specific ORF
#'
#' Gets details of an ORF found within the ORF list based on the start site.
#' This does not perform any checks to see if the ORF is translated, it is simply an in-frame start and stop search.
#' @param orf_object The object created from inspectorORF::find_orfs()
#' @param start_position The position of the start site of interest
#' @return A list containing information about the orfs frame, start and stop codon, and length, or NULL if no ORFs were found
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orfs <- inspectorORF::find_orfs(
#'   tracks,
#'   "ENST00000343702.9"
#' )
#' orf <- inspectorORF::get_orf(
#'   test_orfs,
#'   start_position = 49
#' )
#' orf
#' @export
get_orf <- function(orf_object,
                    start_position)
{
  # Extract the orf(s) matching the requested start position
  temp_result <- lapply(orf_object@orfs, function(x)
  {
    if (!is.null(x$orf_start) && x$orf_start == start_position) x else NULL
  })

  # Discard empty results
  temp_result <- Filter(Negate(is.null), temp_result)

  if (length(temp_result) == 0)
  {
    return(NULL)
  }

  # Return the first matching ORF
  temp_result[[1]]
}

#' Search for ORFs found within the tracks of a transcript
#'
#' Creates a table of ORFs found within a transcript, where ORF is defined as
#' a start (default "AUG", "CUG", "GUG", "UUG", "AAG", and "ACG") and stop ("UAG", "UGA", or "UAA").
#' This does not perform any checks to see if the ORF is translated, it is simply an in-frame start and stop search.
#' @param transcript_tracks The object created from one of the import functions
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param frame_filter The frame(s) to search in (optional) (defailt: 0, 1 and 2)
#' @param start_threshold_filter The minimum number of nucleotides the start site is to be looked at (optional) (default: 17)
#' @param start_codon_filter The start codons to be used for ORF searching (optional) (default: "AUG", "CUG", "GUG", "UUG", "AAG", and "ACG")
#' @param stop_codon_filter The stop codons to be used to determine the end of the ORF (optional) (default: "UAG", "UGA", and "UAA")
#' @param orf_length_filter The minimum length to be considered an ORF (optional) (default: 3)
#' @return An ORF object consisting of frame, start codon, stop codon, start position, stop position, length in nucleotides and length in amino acids
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orfs <- inspectorORF::find_orfs(
#'   tracks,
#'   "ENST00000343702.9"
#' )
#' test_orfs
#' @export
#' @importFrom dplyr case_when mutate filter slice arrange
#' @importFrom tidyr unnest replace_na
find_orfs <- function(transcript_tracks,
                      transcript_filter,
                      frame_filter = c(0, 1, 2),
                      start_threshold_filter = 13,
                      start_codon_filter = c("AUG", "CUG", "GUG", "UUG", "AAG", "ACG"),
                      stop_codon_filter = c("UAG", "UGA", "UAA"),
                      orf_length_filter = 3)
{
  sequence <- transcript_tracks@sequences[[transcript_filter]]
  tracks <- transcript_tracks@tracks[[transcript_filter]] |>
    as.data.frame() |>
    dplyr::filter(name == transcript_tracks@framed_tracks[[1]]) |>
    dplyr::arrange(exon_position)

  # Subset the sequence starting from threshold
  sequence <- substr(sequence, start_threshold_filter, nchar(sequence))

  codon_positions <- seq(1, nchar(sequence) - 2, by = 1)
  codons <- substring(sequence, codon_positions, codon_positions + 2)
  codon_map <- data.frame(
    codon = codons,
    position = codon_positions + start_threshold_filter - 1,
    frame = (codon_positions + start_threshold_filter - 2) %% 3
  )

  codon_map <- codon_map[codon_map$frame %in% frame_filter, ]

  start_codons <- codon_map[codon_map$codon %in% start_codon_filter, ]
  stop_codons  <- codon_map[codon_map$codon %in% stop_codon_filter, ]

  start_codons <- start_codons |>
    dplyr::mutate(is_start = TRUE, name = transcript_filter)

  stop_codons <- stop_codons |>
    dplyr::mutate(is_stop = TRUE, name = transcript_filter)

  tracks_df <- as.data.frame(tracks) |>
    dplyr::filter(name == transcript_tracks@framed_tracks[[1]])

  tracks_df <- tracks_df |>
    dplyr::left_join(start_codons[, c("position", "is_start")], by = c("exon_position" = "position")) |>
    dplyr::left_join(stop_codons[, c("position", "is_stop")], by = c("exon_position" = "position")) |>
    tidyr::replace_na(list(is_start = FALSE, is_stop = FALSE))

  # Get start positions and frames
  start_tracks <- tracks_df |> dplyr::filter(is_start == TRUE)
  start_positions <- start_tracks$exon_position
  start_frames <- start_tracks$framing

  # Get downstream in-frame stops
  stop_positions <- mapply(function(start_pos, frame)
  {
    downstream_stops <- tracks_df |>
      dplyr::filter(exon_position > start_pos, is_stop == TRUE, framing == frame)

    if (nrow(downstream_stops) > 0)
    {
      first_stop <- downstream_stops$exon_position[[1]]

      if ((first_stop - start_pos) >= orf_length_filter)
      {
        return(first_stop)
      }
    }
    return(NA)
  }, start_positions, start_frames, SIMPLIFY = TRUE, USE.NAMES = FALSE)

  valid <- !is.na(stop_positions)
  start_tracks <- start_tracks[valid, ]
  start_positions <- start_positions[valid]
  stop_positions <- stop_positions[valid]

  if (length(stop_positions) == 0)
  {
    return(new("inspectorORF_orflist", orfs = list()))
  }

  orfs <- mapply(
    function(tx_sequence, start_pos, stop_pos, frame, threshold_filter)
    {
      start_codon <- substr(tx_sequence, start_pos - threshold_filter + 1, start_pos - threshold_filter + 3)
      stop_codon  <- substr(tx_sequence, stop_pos- threshold_filter + 1 , stop_pos - threshold_filter + 3)

      .create_orf_list(
        transcript_id = transcript_filter,
        frame = frame,
        start_codon = start_codon,
        stop_codon = stop_codon,
        orf_start = start_pos,
        orf_stop = stop_pos - 1,
        nt_length = stop_pos - start_pos
      )
    },
    sequence,
    start_positions,
    stop_positions,
    start_tracks$framing,
    start_threshold_filter,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  # orfs <- mapply(
  #   .create_orf_list,
  #   transcript_filter,
  #   start_tracks$framing,
  #   substr(sequence, start = start_positions - start_threshold_filter + 1, stop = start_positions - start_threshold_filter + 3),
  #   substr(sequence, start = stop_positions - start_threshold_filter + 1, stop = stop_positions - start_threshold_filter + 3),
  #   start_positions,
  #   stop_positions - 1,
  #   stop_positions - start_positions,
  #   SIMPLIFY = FALSE,
  #   USE.NAMES = FALSE
  # )

  orfs <- orfs[!vapply(orfs, function(x) length(x) == 0 || is.na(x$orf_stop), logical(1))]
  names(orfs) <- paste0(transcript_filter, "_", start_positions, "_", stop_positions - 1)

  new("inspectorORF_orflist", orfs = orfs)
}


# Import reads from a bed file where the tracks are ORF level information
#
# @param bed_file The bed file containing the ORF tracks, each track must contain the trackline of transcript_id_start_stop, e.g ENST00000361764.9_236_400
# @param gtf_file The path to the gtf file
# @param genome_file The path to the genome fasta or 2bit file
# @param framed_tracks The track line consisting of the P-site reads (defaults to p_sites)
#
# @return An inspectorORF_txtracks object consisting of the relevant reads
# @export
#
# @examples
# tracks <- inspectorORF::import_orf_tracks(
#   "test.bed",
#   "gencode.v44.annotation.gtf",
#   "GRCh38.p14.primary_assembly.genome.2bit"
# )
# @importFrom dplyr bind_rows distinct
# @importFrom stringr str_remove
# @importClassesFrom GenomicRanges GRanges
# import_orf_tracks <- function(bed_file,
#                               gtf_file,
#                               genome_file,
#                               framed_tracks = c("p_sites"))
# {
#   bed_tracks <- .import_bed_hack(bed_file)
#
#   names(bed_tracks) <- stringr::str_remove(names(bed_tracks), pattern = "\\_*")
#
#   track_ids <- names(bed_tracks)
#
#   bed_tracks <- do.call(rbind, bed_tracks) |>
#     unique() |>
#     as("GRanges")
#
#   gtf_data <- .import_gtf(gtf_file, track_ids, track_type = "transcript_id")
#
#   transcript_ids <- gtf_data$transcript_id |> unique()
#
#   sequences <- .obtain_sequences(gtf_data, TwoBitFile(genome_file))
#
#   read_names_count <- bed_tracks$name |> unique() |> length()
#
#   tracks <- .get_tracks(bed_tracks, gtf_data, read_names_count, framed_tracks)
#
#   new("inspectorORF_txtracks",
#       transcript_ids = transcript_ids,
#       tracks = tracks,
#       sequences = sequences,
#       framed_tracks = framed_tracks)
# }

#' Used to annotate codons on the ORF plots.
#'
#' A simple wrapper function for a named list to specify codon types you wish to annotate on the ORF/transcript plot
#'
#' @param annotation_codons a vector of strings consisting of which codons to search for within the sequence being plot
#' @param annotate_start a logical indicating if the start codon of the plotted ORF should be anotated
#' @param in_frame a logical indicating if the annotation_codons should be in frame to the ORF or not
#' @param annotate_stop if true and one of annotation_codons or annotate_start is included along with in_frame = TRUE, stop codons will also be annotated. Can be logical value to plot all three stop codons, or vector of strings to specify codons.
#' @param colour a string consisting of the colour to plot the annotated codons in
#'
#' @return a named list of set codon queries for use in the codon_queries option of the orf_plot function
#' @export
#'
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orf_plot <- inspectorORF::orf_plot(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9",
#'   start_position = 456,
#'   codon_queries = list(codon_query(annotate_start = TRUE))
#' )
#' test_orf_plot
codon_query <- function(annotation_codons = NULL, annotate_start = F, in_frame = F, annotate_stop = F, colour = NULL)
{
  if (annotate_start == T && !is.null(annotation_codons))
  {
    stop("Please specify annotate_start and annotation_codons as separate queries.")
  }

  if (annotate_stop == T && in_frame == F && !is.null(annotation_codons))
  {
    stop("annotate_stop requires in_frame to be TRUE when using annotation_codons.")
  }

  if (annotate_stop == T && (annotate_start == F || is.null(annotation_codons)))
  {
    stop("annotate_stop requires either annotate_start or annotation_codons to be set.")
  }

  list(
    annotation_codons = annotation_codons,
    annotate_start = annotate_start,
    in_frame = in_frame,
    annotate_stop = annotate_stop,
    colour = colour
  )
}

#' Plot the reads within an ORF of a transcript
#'
#' @param transcript_tracks The object created from importing the tracks
#' @param orf_object an ORF object returned by inspectorORF::get_orf() - Optional, used instead of transcript_filter, start_position and stop_position
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param start_position The start position to use for the ORF
#' @param stop_position The stop position to use for the ORF (optional)
#' @param plot_region an optional tuple consisting of a region of 5'utr and 3'utr to plot
#' @param plot_colours The custom colour scheme to use for frame 0, 1 and 2 (optional)
#' @param scale_to_psites Should the plot be scaled to the highest P-site peak, therefore cutting off any RNA-Seq reads above this
#' @param plot_transcript_summary Should the remaining P-sites for the full transcript (excluding those in the ORF) be plot
#' @param codon_queries a list consisting of one or more inspectorORF::codon_queries() calls, indiciating any codons to be annotated within the plot
#' @param condition_names Names of any datasets to plot, useful when plotting bed file which consists of reads from multiple datasets
#' @param plot_read_pairs Which RNA-Seq reads are associated with which P-Site reads. See example for further info
#' @param dataset_names Custom naming to display on the plot for any read type
#' @param interactive Enable to return an interactive plotly figure
#' @param one_plot Should a combined figure be returned or individual figures for main plot and any triplet periodicity plots
#' @param legend_position Location of the legend within the figure
#' @param text_size Text size within the figure
#' @return A ggplot of the ORF
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orf_plot <- inspectorORF::orf_plot(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9",
#'   start_position = 456
#' )
#' test_orf_plot
#' @export
orf_plot <- function(
  transcript_tracks,
  orf_object = NULL,
  transcript_filter = NULL,
  start_position = NULL,
  stop_position = NULL,
  plot_region = c(start_position, stop_position),
  plot_colours = c("rna_reads" = "grey60", "0" = "#440854", "1" = "#23A884", "2" = "#FEE725"),
  scale_to_psites = F,
  plot_transcript_summary = F,
  codon_queries = NULL,
  condition_names = c("rna_reads" = ""),
  plot_read_pairs = c("p_sites" = "rna_reads"),
  dataset_names = c("rna_reads" = "RNA-Seq Reads", "p_sites" = "P-Sites"),
  interactive = F,
  one_plot = T,
  legend_position = "bottom",
  text_size = 12
)
{
  .plot_helper(
    transcript_tracks,
    orf_object,
    transcript_filter,
    start_position,
    stop_position,
    plot_region,
    plot_colours,
    scale_to_psites,
    plot_transcript_summary,
    codon_queries,
    condition_names,
    plot_read_pairs,
    dataset_names,
    one_plot,
    interactive,
    legend_position,
    text_size,
    split_exons = F,
    .tx_plot = F
  )
}

#' Extract the nucleotide seqeunce of an ORF within a transcript
#'
#' Extracts the ORF sequence from a transcript track,
#' if the stop position is omitted, the nearest downstream, in-frame stop codon is used from the supplied start codon
#' @param transcript_tracks The object created from inspectorORF::create_tracks()
#' @param orf_object an ORF object returned by inspectorORF::get_orf() - Optional, used instead of transcript_filter, start_position and stop_position
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param start_position The start position to use for the ORF
#' @param stop_position The stop position to use for the ORF (optional)
#' @param stop_codons Which stop codons to consider end of the ORF sequence, defaults to UAG, UAA and UGA
#' @return An RNAStringSet consisting of the ORF nucleotides
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' test_orfs <- inspectorORF::find_orfs(
#'   tracks, "ENST00000343702.9"
#' )
#' inspectorORF::get_orf_nt_seq(
#'   tracks,
#'   test_orfs@orfs[[1]]
#' )
#' @export
#' @importClassesFrom Biostrings RNAStringSet
get_orf_nt_seq <- function(transcript_tracks,
                           orf_object = NULL,
                           transcript_filter = NULL,
                           start_position = NULL,
                           stop_position = NULL,
                           stop_codons = c("UAG", "UAA", "UGA"))
{
  .check_valid_input(orf_object, start_position)

  if (!is.null(orf_object))
  {
    transcript_filter <- orf_object$transcript_id
    start_position <- orf_object$orf_start
    stop_position <- orf_object$orf_stop
  }

  if (!(transcript_filter %in% transcript_tracks@transcript_ids))
  {
    stop(paste("Error!", transcript_filter, "not found within the data"))
  }

  if (is.null(stop_position))
  {
    stop_position <- .get_stop_from_start(transcript_tracks, transcript_filter, start_position, stop_codons)
  }

  if (is.null(stop_position))
  {
    stop("No stop codons found")
  }

  transcript_tracks@sequences[[transcript_filter]] |> substr(start = start_position, stop = stop_position) |>
    Biostrings::RNAStringSet() |>
    setNames(paste(transcript_filter, start_position, stop_position, sep = "_"))
}

#' Extract the amino acid seqeunce of an ORF within a transcript
#'
#' Extracts the ORF sequence from a transcript track,
#' if the stop position is omitted, the nearest downstream, in-frame stop codon is used from the supplied start codon
#' @param transcript_tracks The object created from inspectorORF::create_tracks()
#' @param orf_object an ORF object returned by inspectorORF::get_orf() - Optional, used instead of transcript_filter, start_position and stop_position
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param start_position The start position to use for the ORF
#' @param stop_position The stop position to use for the ORF (optional)
#' @param stop_codons Which stop codons to consider end of the ORF sequence, defaults to UAG, UAA and UGA
#' @return An AAStringSet consisting of the ORF amino acid residues
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' orf_aa_sequence <- inspectorORF::get_orf_aa_seq(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9",
#'   start_position = 456
#' )
#' orf_aa_sequence
#' @export
#' @importFrom bioseq seq_translate as_rna
#' @importClassesFrom Biostrings AAStringSet
get_orf_aa_seq <- function(transcript_tracks,
                           orf_object = NULL,
                           transcript_filter,
                           start_position = NULL,
                           stop_position = NULL,
                           stop_codons = c("UAG", "UAA", "UGA"))
{
  bioseq::seq_translate(
    bioseq::as_rna(
      as.vector(
        get_orf_nt_seq(
          transcript_tracks,
          orf_object,
          transcript_filter,
          start_position,
          stop_position,
          stop_codons
        )
      )
    )
  ) |> Biostrings::AAStringSet()
}

#' Obtain reads for all three frames of an ORF
#'
#' Obtain the RNA-Seq and Ribo-Seq reads for all three frames from the specified ORF
#' @param transcript_tracks the inspectorORF tracks object generated by import_transcript_bed or get_transcript_tracks
#' @param orf_object an ORF object returned by inspectorORF::get_orf() - Optional, used instead of transcript_filter, start_position and stop_position
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param start_position The start position to use for the ORF
#' @param stop_position The stop position to use for the ORF (optional)
#' @param stop_codons Which stop codons to consider end of the ORF sequence, defaults to UAG, UAA and UGA
#' @return a table consisting of summarised reads for each frame
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' orf_framing <- inspectorORF::get_orf_framing(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9",
#'   start_position = 456
#' )
#' orf_framing
#' @export
#' @importFrom dplyr filter summarise group_by
get_orf_framing <- function(transcript_tracks,
                            orf_object = NULL,
                            transcript_filter = NULL,
                            start_position = NULL,
                            stop_position = NULL,
                            stop_codons = c("UAG", "UAA", "UGA"))
{
  # .check_valid_input(orf_object, start_position)

  if (!is.null(orf_object))
  {
    transcript_filter <- orf_object$transcript_id
    start_position <- orf_object$orf_start
    stop_position <- orf_object$orf_stop
  }

  if (!(transcript_filter %in% transcript_tracks@transcript_ids))
  {
    stop(paste("Error!", transcript_filter, "not found within the data"))
  }

  if (is.null(stop_position))
  {
    stop_position <- .get_stop_from_start(transcript_tracks, transcript_filter, start_position, stop_codons)
  }

  transcript_tracks@tracks[[transcript_filter]] |>
    dplyr::filter(name == transcript_tracks@framed_tracks &
                    exon_position > (start_position - 1) &
                    exon_position < (stop_position + 1)) |>
    dplyr::group_by(framing) |>
    dplyr::summarise(totals = sum(score))
}
