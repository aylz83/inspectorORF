#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom rtracklayer TwoBitFile
#'
#' @slot tracks A GRanges object with track data.
#' @slot track_ids A character vector of track IDs (typically either gene_ids or transcript_ids).
#' @slot gtf A GRanges object representing a gene annotation (e.g., GTF).
#' @slot genome_file A TwoBitFile representing the genome sequence.
#' @slot framed_tracks A character vector of tracks marked as framed.
setClass(
  "inspectorORF_gtracks",
  slots = c(
    tracks = "GRanges",
    track_ids = "vector",
    gtf = "GRanges",
    genome_file = "TwoBitFile",
    framed_tracks = "vector",
    extra_data = "list"
  ),
  prototype = list(
    tracks = GRanges(),
    track_ids = c(),
    gtf = GRanges(),
    genome_file = NULL,
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
#' tracks <- inspectorORF::import_gene_tracks("test.bed", "gencode.v44.nnotation.gtf", "GRCh38.p14.primary_assembly.genome.2bit")
#' @importFrom dplyr bind_rows distinct
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom rtracklayer TwoBitFile
#' @importClassesFrom rtracklayer TwoBitFile
import_gene_tracks <- function(bed_file, gtf_file, genome_file, framed_tracks = c("p_sites"))
{
  bed_tracks <- .import_bed_hack(bed_file)

  track_ids <- names(bed_tracks)

  bed_tracks <- dplyr::bind_rows(bed_tracks) |> dplyr::distinct() |> as("GRanges")

  gtf_data <- .import_gtf(gtf_file, track_ids)

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
    genome_file = TwoBitFile(genome_file),
    # read_names = read_names,
    framed_tracks = framed_tracks,
    # extra_data = extra_data
  )
}
