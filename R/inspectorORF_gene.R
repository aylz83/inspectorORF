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

