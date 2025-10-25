#' @importFrom GenomicRanges GRangesList
#' @importClassesFrom GenomicRanges GRangesList
#'
#' @slot transcript_ids A character vector of transcript IDs.
#' @slot tracks A GRangesList of genomic tracks per transcript.
#' @slot sequences A character vector of associated sequences.
#' @slot framed_tracks A character vector of framed track IDs.
setClass(
  "inspectorORF_txtracks",
  slots = c(
    transcript_ids = "vector",
    tracks = "GRangesList",
    sequences = "vector",
    framed_tracks = "vector"
  ),
  prototype = list(
    transcript_ids = c(),
    tracks = GRangesList(),
    sequences = c(),
    framed_tracks = c()
  )
)

# .add_introns <- function(tracks, read_names)
# {
#   introns_to_add <- tracks@tracks |> filter(at_exon_end == T) |> pull(introns_to_add, exon_position)
#
#   if (length(introns_to_add) == 1)
#   {
#     names(introns_to_add) <- as.numeric(names(introns_to_add)) + introns_to_add
#   }
#   else
#   {
#     for (at in 2:length(introns_to_add))
#     {
#       names(introns_to_add)[at] <- as.numeric(names(introns_to_add)[at]) + sum(introns_to_add[1:at - 1])
#     }
#   }
#
#
#   final_data <- tracks@tracks
#
#   for (position in seq_along(introns_to_add))
#   {
#     # print(paste("Adding", introns_to_add[position], "introns at position", names(introns_to_add)[position]))
#     final_data <- final_data |> add_row(is_exon = rep(F, introns_to_add[position]),
#                                          .after = as.numeric(names(introns_to_add)[position]))
#   }
#
#   final_data |> mutate(with_introns_position = row_number()) |>
#     mutate_at(tracks@read_names, ~replace_na(., 0))
# }

# .splice_plot <- function(transcript_data)
# {
#   # frame <- transcript_data |> dplyr::filter(exon_position == start_position) |> dplyr::pull(framing)
#
#   intron_junctions <- transcript_data |> filter(at_exon_end == T)
#   intron_count <- intron_junctions |> nrow()
#
#   exon_junctions <- lapply(1:intron_count, function(intron_junction)
#   {
#     intron_position <- intron_junctions |> dplyr::slice(intron_junction) |> pull(with_introns_position)
#     exon_junction <- transcript_data |> dplyr::slice(-c(1:intron_position)) |>
#       dplyr::filter(is_exon == T) |>
#       dplyr::slice(1)
#   }) |> bind_rows(intron_junctions,
#                    transcript_data |> dplyr::slice(1),
#                    transcript_data |> dplyr::slice(n())) |>
#     dplyr::arrange(with_introns_position) |>
#     dplyr::mutate(with_introns_position = as.factor(as.character(with_introns_position))) |>
#     dplyr::mutate(exon_feature = as.factor(rep(c("start", "end"), length.out = n()))) |>
#     dplyr::select(exon_number, with_introns_position, exon_feature) |>
#     pivot_wider(names_from = exon_feature, values_from = with_introns_position, values_fn = list) |>
#     unnest(cols = c(start, end))
#
#   exon_junctions |>
#     ggplot() +
#     geom_hline(yintercept = 0.5) +
#     geom_rect(aes(xmin = start, xmax = end,
#                   ymin = 0, ymax = 1),
#               color = "black") +
#     scale_x_discrete(limits=factor(1:nrow(transcript_data)), labels = as.character(1:nrow(transcript_data))) +
#     theme_void() +
#     theme(legend.position = "none")
# }

#' Get a vector of all transcript ids found within a particular gene_id
#'
#' @param gene_tracks The gene tracks object, obtained from import_gene_tracks
#' @param gene_filter The gene_id of interest
#'
#' @return A vector of transcript_ids associated with the gene_id
#' @export
#'
#' @examples
#' gene_tracks <- inspectorORF::import_gene_bed(
#'   system.file("example_data", "example_gene.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF")
#' )
#' print(inspectorORF::get_transcripts_for_gene(
#'   gene_tracks,
#'   gene_filter = "ENSG00000074527.1"
#' ))
#' @importFrom dplyr filter pull
get_transcripts_for_gene <- function(gene_tracks, gene_filter)
{
  gene_tracks@gtf |> as.data.frame() |>
    dplyr::filter(gene_id == gene_filter) |>
    dplyr::pull(transcript_id) |>
    unname()
}

#' Import reads from a bed file where the tracks are transcript level information
#'
#' @param bed_file The bed file containing the ORF tracks, each track must contain the trackline of a transcript_id
#' @param gtf_file The path to the gtf file
#' @param genome_file The path to the genome fasta or 2bit file
#' @param framed_tracks The track line consisting of the P-site reads (defaults to p_sites)
#'
#' @return An inspectorORF_txtracks object consisting of the relevant reads
#' @export
#'
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF")
#' )
#' @importFrom dplyr bind_rows distinct
#' @importClassesFrom rtracklayer TwoBitFile
#' @importFrom methods as new
import_transcript_bed <- function(bed_file,
                                  gtf_file,
                                  genome_file,
                                  framed_tracks = c("p_sites"))
{
  bed_tracks <- .import_bed_hack(bed_file)

  track_ids <- names(bed_tracks)

  bed_tracks <- dplyr::bind_rows(bed_tracks) |> dplyr::distinct() |> as("GRanges")

  gtf_data <- .import_gtf(gtf_file, track_ids, track_type = "transcript_id")

  transcript_ids <- gtf_data$transcript_id |> unique()

  sequences <- .obtain_sequences(gtf_data, rtracklayer::TwoBitFile(genome_file))

  read_names_count <- bed_tracks$name |> unique() |> length()

  tracks <- .get_tracks(bed_tracks, gtf_data, read_names_count, framed_tracks)

  new("inspectorORF_txtracks",
      transcript_ids = transcript_ids,
      tracks = tracks,
      sequences = sequences,
      framed_tracks = framed_tracks)
}

#' Extracts reads from gene level data relevant to a transcript of interest
#'
#' @param gene_tracks The gene tracks object obtained from inspectorORF::import_gene_tracks
#' @param transcript_filter The transcript ids of interest
#'
#' @return A transcript tracks object
#' @export
#'
#' @examples
#' gene_tracks <- inspectorORF::import_gene_bed(
#'   system.file("example_data", "example_gene.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF")
#' )
#' tx_tracks <- inspectorORF::gene_to_transcript_tracks(
#'   gene_tracks,
#'   transcript_filter = c("ENST00000343702.9", "ENST00000344911.8")
#' )
#' @importFrom plyranges filter
gene_to_transcript_tracks <- function(gene_tracks,
									  transcript_filter)
{
  gtf_subset <- gene_tracks@gtf |>
    plyranges::filter(transcript_id %in% transcript_filter)

  transcript_ids <- gtf_subset$transcript_id |> unique()

  sequences <- .obtain_sequences(gtf_subset, gene_tracks@genome_file)

  read_names_count <- gene_tracks@tracks$name |> unique() |> length()

  tracks <- .get_tracks(gene_tracks@tracks, gtf_subset, read_names_count, gene_tracks@framed_tracks)

  new("inspectorORF_txtracks",
      transcript_ids = transcript_ids,
      tracks = tracks,
      sequences = sequences,
      framed_tracks = gene_tracks@framed_tracks)
}

#' Plot the reads for a full transcript
#'
#' Plots the Ribo-Seq, RNA-Seq, and Mass-Spec data across a full transcript.
#' Optionally, with an ORF of interest outlined.
#' @param transcript_tracks The object created from inspectorORF::create_tracks()
#' @param orf_object an ORF object returned by inspectorORF::get_orf() - Optional, used instead of transcript_filter, start_position and stop_position
#' @param transcript_filter A character vector consisting of your transcript ID of interest
#' @param start_position The start location of the ORF to highlight (optional)
#' @param stop_position The stop location of the ORF to highlight, if omitted, the nearest downstream, in-frame stop codon from the supplied start position will be used.
#' @param plot_colours The colour scheme for frame 0, 1 and 2 (optional).
#' @param scale_to_psites Should the plot be scaled to the highest P-site peak, therefore cutting off any RNA-Seq reads above this
#' @param plot_transcript_summary Should the remaining P-sites for the full transcript (excluding those in the ORF) be plot
#' @param codon_queries an optional list consisting of one or more inspectorORF::codon_queries() calls, indiciating any codons to be annotated within the plot
#' @param condition_names Names of any datasets to plot, useful when plotting bed file which consists of reads from multiple datasets
#' @param plot_read_pairs Which RNA-Seq reads are associated with which P-Site reads. See example for further info
#' @param dataset_names Custom naming to display on the plot for any read type
#' @param interactive Enable to return an interactive plotly figure
#' @param one_plot Should a combined figure be returned or individual figures for main plot and any triplet periodicity plots
#' @param legend_position Location of the legend within the figure
#' @param text_size Text size within the figure
#' @return a ggplot object.
#' @export
#'
#' @examples
#' tracks <- inspectorORF::import_transcript_bed(
#'   system.file("example_data", "example_tracks.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF"),
#' )
#' tx_plot <- inspectorORF::transcript_plot(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9"
#' )
#'
#' # Picking out an ORF of interest within the transcript plot
#' tx_plot <- inspectorORF::transcript_plot(
#'   tracks,
#'   transcript_filter = "ENST00000343702.9",
#'   start_position = 10,
#'   stop_position = 30
#' )
transcript_plot <- function(transcript_tracks,
                            orf_object = NULL,
                            transcript_filter = NULL,
                            start_position = NULL,
                            stop_position = NULL,
                            plot_colours = c("rna_reads" = "grey60", "0" = "#440854", "1" = "#23A884", "2" = "#FEE725"),
                            scale_to_psites = F,
                            plot_transcript_summary = F,
                            codon_queries = NULL,
                            condition_names = c("rna_reads" = ""),
                            plot_read_pairs = c("p_sites" = "rna_reads"),
                            dataset_names = c("rna_reads" = "RNA-Seq Reads",
                                              "p_sites" = "P-Sites"),
                            one_plot = T,
                            interactive = F,
                            legend_position = "bottom",
                            text_size = 12)
{
  .plot_helper(transcript_tracks,
           orf_object,
           transcript_filter,
           start_position,
           stop_position,
           plot_region = c(-1, -1),
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
           .tx_plot = T)
}
