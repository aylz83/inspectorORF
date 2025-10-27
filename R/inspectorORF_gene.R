#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importClassesFrom rtracklayer TwoBitFile
#'
#' @slot tracks A GRanges object with track data.
#' @slot track_ids A character vector of gene_ids.
#' @slot gtf A GRanges object representing a gene annotation (e.g., GTF).
#' @slot genome_file A TwoBitFile representing the genome sequence.
#' @slot framed_tracks A character vector of tracks marked as framed.
setClass(
  "inspectorORF_gtracks",
  slots = c(
    tracks = "GRangesList",
    track_ids = "vector",
    gtf = "GRanges",
    genome_file = "TwoBitFile",
    framed_tracks = "vector"
  ),
  prototype = list(
    tracks = GRangesList(),
    track_ids = c(),
    gtf = GRanges(),
    genome_file = new("TwoBitFile"),
    framed_tracks = c()
  )
)

#' Import reads from a bed file where the tracks are gene level information
#'
#' @param bed_file The bed file containing the gene tracks, each track must contain the trackline of gene_id
#' @param gtf_file The path to the gtf file
#' @param genome_file The path to the genome fasta or 2bit file
#' @param framed_tracks The track line consisting of the P-site reads (defaults to p_sites)
#'
#' @return An inspectorORF_gtracks object consisting of the relevant reads
#' @export
#'
#' @examples
#' gene_tracks <- inspectorORF::import_gene_bed(
#'   system.file("example_data", "example_gene.bed", package = "inspectorORF"),
#'   gtf_file = system.file("example_data", "annotation_subset.gtf", package = "inspectorORF"),
#'   genome_file = system.file("example_data", "chr12.2bit", package = "inspectorORF")
#' )
#' @importFrom dplyr bind_rows distinct
#' @importFrom methods as
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom rtracklayer TwoBitFile
import_gene_bed <- function(bed_file, gtf_file, genome_file, framed_tracks = c("p_sites"))
{
  bed_tracks <- .import_bed_hack(bed_file, as_ranges = TRUE)

  track_ids <- names(bed_tracks)

  gtf_data <- .import_gtf(gtf_file, track_ids, track_type = "gene_id")

  # if (is.null(extra_data_file))
  # {
  #   extra_data <- list()
  # }
  # else
  # {
  #   extra_data <- .import_bed_hack(extra_data_file)
  # }

  new(
    "inspectorORF_gtracks",
    tracks = bed_tracks,
    track_ids = track_ids,
    gtf = gtf_data,
    genome_file = rtracklayer::TwoBitFile(genome_file),
    # read_names = read_names,
    framed_tracks = framed_tracks
    # extra_data = extra_data
  )
}

.gene_plot <- function(
  gene_tracks,
  gene_filter,
  transcript_filter = NULL,
  orf_object = NULL,
  orf_start = NULL,
  orf_end = NULL,
  plot_colours = c("rna_reads" = "grey60", "0" = "#440854", "1" = "#23A884", "2" = "#FEE725"),
  scale_to_psites = F,
  plot_transcript_summary = F,
  shrink_introns = F,
  codon_queries = NULL,
  condition_names = c("rna_reads" = ""),
  plot_read_pairs = c("p_sites" = "rna_reads"),
  dataset_names = c("rna_reads" = "RNA-Seq Reads", "p_sites" = "P-Sites"),
  one_plot = T,
  interactive = F,
  legend_position = "bottom",
  text_size = 12
)
{
	track_to_plot <- gene_tracks@tracks[[gene_filter]] |>
		as.data.frame() |>
		dplyr::distinct(seqnames, start, end, strand, name, score, .keep_all = TRUE)

	read_names_count <- track_to_plot$name |> unique() |> length()

	track_to_plot <- track_to_plot |>
		dplyr::mutate(framing = as.factor(ifelse(name %in% gene_tracks@framed_tracks, rep(c(0, 1, 2), each = read_names_count, length.out = dplyr::n()), name)))

	p_sites <- track_to_plot |> dplyr::select(seqnames, start, end, strand, name, score, framing) |>
		dplyr::filter(name %in% names(plot_read_pairs)) |>
		dplyr::rename(p_sites = score,
					  p_site_framing = framing) |>
		dplyr::arrange(match(name, names(plot_read_pairs))) |>
		dplyr::mutate(name = rep(plot_read_pairs, each = length(name) / length(plot_read_pairs)))

	# print(track_to_plot)
	track_to_plot <- track_to_plot |> dplyr::left_join(p_sites, by = c("seqnames", "start", "end", "strand", "name")) |>
		dplyr::filter(!(name %in% names(plot_read_pairs)))

	filtered_framed_tracks <- intersect(gene_tracks@framed_tracks, names(plot_read_pairs))

	plot_labels <- dataset_names[!(names(dataset_names) %in% filtered_framed_tracks)]
	plot_labels[names(plot_labels)] <- ifelse(
		names(plot_labels) %in% names(condition_names) & condition_names[names(plot_labels)] != "",
		paste0(condition_names[names(plot_labels)], "\n", dataset_names[names(plot_labels)]),
		paste0(dataset_names[names(plot_labels)])
	)

	plot_labels[plot_read_pairs] <- paste0(plot_labels[plot_read_pairs], "\n+", dataset_names[names(plot_read_pairs)])

	legend_labels <- sapply(names(plot_colours), function(nm)
	{
		dataset_label <- if (!is.null(dataset_names) && nm %in% names(dataset_names))
		{
			dataset_names[[nm]]
		}
		else
		{
			nm
		}

		label <- if (!is.null(condition_names) && nm %in% names(condition_names) && condition_names[[nm]] != "")
		{
			paste0(condition_names[[nm]], "\n", dataset_label)
		}
		else
		{
			dataset_label
		}

		if (nm %in% c("0", "1", "2"))
		{
			paste("Frame", nm)
		}
		else
		{
			label
		}
	}, USE.NAMES = TRUE)

	plot_result <- ggplot(
		track_to_plot,
		aes(y = score,
			x = start,
			colour = framing,
			fill = framing)
	) +
		geom_bar(stat = "identity", na.rm = TRUE) +
		geom_bar(
			aes(x = start,
				y = p_sites,
				colour = p_site_framing,
				fill = p_site_framing),
			stat = "identity",
			na.rm = TRUE
		) +
		scale_color_manual(
			values = plot_colours,
			aesthetics = c("fill", "colour"),
			labels = legend_labels
		) +
		scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
		coord_cartesian(expand = FALSE) +
		xlab("Position") +
		ylab(NULL) +
		theme_minimal() +
		theme(
			legend.position = legend_position,
			legend.title = element_blank(),
			panel.spacing = grid::unit(10, "pt"),
			# panel.border = element_rect(colour = "black", fill = NA, linewidth = .5),
			panel.background = element_blank(),
			# axis.line = element_line(colour = "black"),
			# axis.title.y = element_text(color = "grey30"),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank(),
			strip.background.y = element_blank(),
			strip.placement = "left",
			# strip.text = element_text(size = 12),
			strip.text.y.left = element_text(angle = 0),
			axis.ticks.length = grid::unit(0, "points"),
			text = element_text(size = text_size),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)
}
