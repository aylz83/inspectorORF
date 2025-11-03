.check_files_exist <- function(files, label)
{
  if (any(grepl("\\.bed\\.gz$", files)) && !requireNamespace("R.utils", quietly = TRUE))
  {
    stop("R.utils is required for opening gzip compressed bed files.")
  }

  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0)
  {
    stop(paste0("The following ", label, " files do not exist:\n",
      paste(missing_files, collapse = "\n")
    ))
  }
}

.debug_inform <- function(obj, output_file = "debug_log.txt")
{
  # output <- capture.output(print(obj))
  # msg <- paste(output, collapse = "\n")
  # log_entry <- paste0(format(Sys.time(), "%H:%M:%S"), ": ", msg, "\n")
  # cat(log_entry, file = output_file, append = TRUE)
}

.read_big_text <- function(x)
{
  f = file(x, "rb")
  a = readChar(f, file.info(x)$size, useBytes = T)
  a <- strsplit(a, "\n", fixed = T ,useBytes = T)[[1]]
  close(f)
  return(a)
}

.extract_attributes <- function(gtf_attributes, att_of_interest)
{
  att <- unlist(strsplit(gtf_attributes, " "))
  gsub("\"|;","", att[which(att %in% att_of_interest) + 1])
}

#' @importFrom plyranges filter
.filter_gtf_df <- function(gtf_data, ids, attribute_column, track_type)
{
  # if (track_type == "gene_id")
  # {
  #   filtered <- paste(ids, collapse = "|")
  #   return (gtf_data |> plyranges::filter(type == "exon" & grepl(filtered, !!as.name(attribute_column))))
  # }

  if (length(ids) == 0 || is.null(ids) || is.na(ids))
  {
    return(gtf_data)
  }

  filtered <- paste(ids, collapse = "|")
  gtf_data |> plyranges::filter(type == "exon" & grepl(filtered, !!as.name(attribute_column)))
}

#' @importFrom plyranges select
.select_gtf_df <- function(gtf_data)
{
  gtf_data |> plyranges::select(
    gene_id,
    gene_name,
    chr,
    start,
    end,
    strand,
    transcript_id,
    transcript_name,
    transcript_type,
    exon_number
  )
}

#' @importFrom plyranges mutate
.process_gtf_table <- function(gtf_data)
{
  gtf_data |> plyranges::mutate(
    gene_id = sapply(attributes, .extract_attributes, "gene_id"),
    gene_name = sapply(attributes, .extract_attributes, "gene_name"),
    transcript_id = sapply(attributes, .extract_attributes, "transcript_id"),
    transcript_name = sapply(attributes, .extract_attributes, "transcript_name"),
    transcript_type = sapply(attributes, .extract_attributes, "transcript_type"),
    exon_number = as.integer(sapply(attributes, .extract_attributes, "exon_number"))
  )
}

# I've found importing GTF can be faster if i process this myself due to the
# ability of fread to use multiple cores
#' @importFrom data.table fread
#' @importClassesFrom GenomicRanges GRanges
.import_gtf_hack <- function(gtf_file, filter, track_type)
{
  data.table::fread(
    gtf_file,
    sep = "\t",
    header = F,
    col.names = c(
      "chr",
      "source",
      "type",
      "start",
      "end",
      "score",
      "strand",
      "phase",
      "attributes"
    )
  ) |> .filter_gtf_df(filter, "attributes", track_type) |>
    .process_gtf_table() |>
    .select_gtf_df() |>
    as("GRanges")
}

#' @importFrom data.table fread
.check_id_type <- function(gtf_file, filter)
{
  rows <- data.table::fread(
    gtf_file,
    sep = "\t",
    header = F,
    nrows = 10,
    col.names = c(
      "chr",
      "source",
      "type",
      "start",
      "end",
      "score",
      "strand",
      "phase",
      "attributes"
    )
  ) |> .filter_gtf_df(filter, "attributes", track_type) |>
    .process_gtf_table()

  if (filter %in% rows$gene_id)
  {
    return("gene_id")
  }
  else if (filter %in% rows$transcript_id)
  {
    return("transcript_id")
  }
  else
  {
    return("unknown")
  }
}

#' @importClassesFrom GenomicRanges GRanges
.import_gtf <- function(gtf_file, track_ids, track_type = "gene_id")
{
  if (c("GRanges") %in% class(gtf_file) == T)
  {
    gtf_data <- gtf_file

    gtf_data <- gtf_data |> .filter_gtf_df(track_ids, track_type)
    gtf_data <- gtf_data |> .select_gtf_df()
  }
  else if (any(class(gtf_file) %in% c("data.frame", "data.table")) == T)
  {
    if ("attributes" %in% colnames(gtf_file))
    {
      filter_column <- "attributes"
    }
    else
    {
      filter_column <- track_type
    }

    gtf_data <- gtf_file |>
      .filter_gtf_df(track_ids, filter_column)

    if ("attributes" %in% colnames(gtf_file))
    {
      gtf_data <- gtf_data |> .process_gtf_table()
    }

    gtf_data <- gtf_data |> .select_gtf_df()
  }
  else
  {
    gtf_data <- .import_gtf_hack(gtf_file, track_ids, track_type)
  }

  gtf_data
}

#' @importFrom data.table fread
.import_bed_hack <- function(bed_file, as_ranges = FALSE)
{
  big_data <- .read_big_text(bed_file)
  track_positions <- which(grepl("^track name=", big_data))
  track_names <- sub("^track name=", "", big_data[track_positions])

  track_pairs <- sort(c(track_positions[-c(1)] - 1, (track_positions + 1), length(big_data)))

  results <- setNames(lapply(seq(from = 1, to = length(track_pairs / 2), by = 2), function(at)
  {
    data <- data.table::fread(text = big_data[track_pairs[at]:track_pairs[at + 1]],
          col.names = c("seqnames", "start", "end", "name", "score", "strand")) |>
      dplyr::mutate(start = start + 1)# |>

    if (as_ranges == TRUE)
    {
      data |> as("GRanges")
    }
    else
    {
      data
    }
  }), track_names)

  if (as_ranges)
  {
    results |> as("GRangesList")
  }
  else
  {
    results
  }
}

#' @importFrom dplyr mutate distinct select rename rename_with group_by ungroup
#' @importFrom data.table fread
.import_coverage_bed <- function(bed_file, score_name, n_check = 10)
{
  # Read first few lines to check structure
  sample_data <- data.table::fread(bed_file, nrows = n_check, header = FALSE)
  n_cols <- ncol(sample_data)

  if (n_cols >= 8)
  {
    # Full BED with extra columns generated by bedtools coverage
    col_names <- c("seqnames", "start", "end", "biotype", "gene_id", "strand", "position", "score")
    bed_data <- data.table::fread(bed_file, col.names = col_names) |>
      dplyr::mutate(chr_pos = start + position) |>
      dplyr::distinct(seqnames, chr_pos, strand, .keep_all = TRUE) |>
      dplyr::select(seqnames, chr_pos, strand, score) |>
      dplyr::rename_with(~ score_name, .cols = "score") |>
      dplyr::rename(start = chr_pos)

  }
  else if (n_cols >= 5)
  {
    # Minimal BED: seqnames, start, end, strand, score
    col_names <- c("seqnames", "start", "end", "strand", "score")
    bed_data <- data.table::fread(bed_file, col.names = col_names)

    # Check if 'start' is already unique/monotonic within seqnames
    bed_sample_check <- bed_data[1:min(n_check, nrow(bed_data)), ]

    needs_increment <- any(
      duplicated(bed_sample_check$start)
    )

    if (needs_increment)
    {
      # Manual increment: add 0:(nrow-1) within each seqname
      bed_data <- bed_data |>
        dplyr::group_by(seqnames) |>
        dplyr::mutate(start = start + dplyr::row_number() - 1) |>
        dplyr::ungroup()
    }

    # Rename score
    bed_data <- bed_data |> dplyr::rename_with(~ score_name, .cols = "score")

  }
  else
  {
    stop("BED file has insufficient columns. Expected at least 5 columns.")
  }

  return(bed_data)
}

# from:
# https://github.com/mdeber/BRGenomics/blob/dbab7113a1556c82053880cc0f106e97ac397f34/R/dataset_functions.R#L106
# borrowed from BRGenomics (https://github.com/mdeber/BRGenomics) due to this package not installing via BiocManager anymore
# A copy of the licence can be found in LICENSE.BRGenomics and
#' @importFrom methods is
#' @importFrom GenomicRanges width GPos mcols mcols<- isDisjoint findOverlaps
#'   GRanges sort
#' @importFrom S4Vectors from
#' @importFrom parallel mclapply
.brgenomics_makeGRangesBRG <- function(dataset.gr, ncores = getOption("mc.cores", 2L))
{
  if (is.list(dataset.gr) || is(dataset.gr, "GRangesList"))
    return(parallel::mclapply(dataset.gr, .brgenomics_makeGRangesBRG, mc.cores = ncores))

  if (!GenomicRanges::isDisjoint(dataset.gr))
    stop("Input dataset.gr is not disjoint. See documentation")

  # separate wide granges
  is_wide <- GenomicRanges::width(dataset.gr) > 1L
  gr_wide <- dataset.gr[is_wide]

  # make width 1, reverse map, and add metadata
  gp <- GenomicRanges::GPos(gr_wide)
  hits <- GenomicRanges::findOverlaps(gr_wide, gp)
  GenomicRanges::mcols(gp) <- GenomicRanges::mcols(gr_wide)[S4Vectors::from(hits), , drop = FALSE]

  # combine and sort
  GenomicRanges::sort(c(dataset.gr[!is_wide], GRanges(gp)))
}

#' @importFrom dplyr mutate select rename
.import_orfquant <- function(for_file, score_name, orfquant_results_type)
{
  loaded_name <- load(file = for_file)
  get(loaded_name)[[orfquant_results_type]] |> .brgenomics_makeGRangesBRG() |>
    as.data.frame() |>
    dplyr::select(seqnames, start, strand, score) |>
    dplyr::mutate(start = as.numeric(start),
                  start = as.numeric(start)) |>
    dplyr::rename_with(~ score_name, .cols = "score")
}

#' @importFrom dplyr mutate arrange group_by summarise
#' @importFrom tibble deframe
#' @importClassesFrom Biostrings RNAStringSet
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom rtracklayer getSeq
.obtain_sequences <- function(gtf_data, genome_file)
{
  gtf_data |> as.data.frame() |>
    dplyr::mutate(
      sequence = rtracklayer::getSeq(
        genome_file,
        which = GenomicRanges::GRanges(paste0(seqnames, ":", start, "-", end, ":", strand))) |>
      Biostrings::RNAStringSet() |>
      as.character()
    ) |>
    dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), xtfrm(start))) |>
    dplyr::group_by(transcript_id) |>
    dplyr::summarise(sequence = paste(sequence, collapse = "")) |>
    tibble::deframe()
}

#' @importFrom dplyr group_by mutate ungroup n
#' @importFrom methods as
#' @importFrom plyranges join_overlap_inner_within_directed
#' @importClassesFrom GenomicRanges GRanges
.get_tracks <- function(tracks, gtf_subset, read_names_count, framed_tracks)
{
  tracks <- plyranges::join_overlap_inner_within_directed(tracks, gtf_subset, maxgap = -1L, minoverlap = 0L) |>
    as.data.frame() |> # faster to convert to and from dataframe than it is to use plyranges::group_by?!
    dplyr::distinct(seqnames, start, end, strand, transcript_id, name, score, .keep_all = TRUE) |>
    dplyr::group_by(transcript_id) |>
    dplyr::mutate(
      is_intron = grepl("intron", feature_number),

      exon_row = ifelse(
        is_intron, NA_integer_,
        cumsum(!is_intron & !duplicated(paste(transcript_id, feature_number, start)))
      ),

      exon_frame = ifelse(is_intron, NA_integer_, (exon_row - 1) %% 3),

      genomic_position = cumsum(!duplicated(paste(transcript_id, start, feature_number))),

      og_framing = dplyr::case_when(
        is_intron & name %in% framed_tracks ~ "intron_psite",
        is_intron ~ "intron",
        TRUE ~ as.character(exon_frame)
      ),

      framing = dplyr::case_when(
      	is_intron & name %in% framed_tracks ~ "intron_psite",
      	is_intron ~ "intron",
        name %in% framed_tracks ~ as.character(exon_frame),
        TRUE ~ name
      )
      # og_framing = as.factor(rep(c(0, 1, 2), each = read_names_count, length.out = dplyr::n())),
      # framing = as.factor(ifelse(name %in% framed_tracks, rep(c(0, 1, 2), each = read_names_count, length.out = dplyr::n()), name))
      # at_exon_end = ifelse(strand == "+", start == end_exon & row_number() > read_names_count, start == start_exon & row_number() < (length(new_tracks) - read_names_count)),
      # introns_to_add = ifelse(at_exon_end == F, NA, round(abs(start_exon - end_exon) / 10)),
      # is_exon = T
    ) |>
    dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), xtfrm(start))) |>
    dplyr::ungroup() |>
    dplyr::select(
      seqnames,
      start,
      end,
      strand,
      transcript_id,
      name,
      score,
      is_intron,
      genomic_position,
      og_framing,
      framing,
      feature_number
    ) |>
    as("GRanges")

  split(tracks, tracks$transcript_id) |> as("GRangesList")
}

#' @importFrom stats setNames
.cat_lists <- function(list1, list2)
{
  keys <- unique(c(names(list1), names(list2)))
  setNames(mapply(c, list1[keys], list2[keys], SIMPLIFY = FALSE), keys)
}

.get_stop_from_start <- function(inspectorORF_tracks, transcript_filter, start_position, stop_codons)
{
  sequence <- inspectorORF_tracks@sequences[[transcript_filter]]

  track_df <- as.data.frame(inspectorORF_tracks@tracks[[transcript_filter]])
  track_df <- track_df[track_df$name == inspectorORF_tracks@framed_tracks[[1]], ]

  # obtain the frame for the specified start position
  orf_frame <- as.numeric(track_df$framing[start_position])

  # filter tracks to remove everything before the start position of the ORF
  remaining <- track_df[(start_position):nrow(track_df), ]
  remaining <- remaining[as.numeric(remaining$framing) == orf_frame, ]

  # loop through positions in-frame and look for stop codon
  for (i in seq_len(nrow(remaining)))
  {
    codon_start <- remaining$genomic_position[i]
    codon <- substr(sequence, codon_start, codon_start + 2)

    if (codon %in% stop_codons)
    {
      return(codon_start - 1)
    }
  }

  # no stop codon found
  return(NULL)
}

#' @importFrom ggplot2 ggplot geom_bar aes scale_color_manual coord_cartesian xlab ylab theme_minimal theme facet_grid geom_blank scale_y_continuous labeller
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggh4x facetted_pos_scales
#' @importFrom grid unit arrow
.codon_plot <- function(
  plot_result,
  codon_query,
  row_sizes
)
{
  # Compute minimum spacing between sorted exon positions
  codon_positions <- sort(codon_query$genomic_position)
  spacing <- diff(codon_positions)

  # Count how many pairs are "close"
  close_threshold <- 5
  num_close_pairs <- sum(spacing < close_threshold)

  # Scale panel height by number of close pairs
  min_codon_row_size <- 0.1
  max_codon_row_size <- 1
  codon_row_size <- min(
    max(
      min_codon_row_size,
      0.1 + 0.08 * num_close_pairs  # adjust scaling factor as needed
    ),
    max_codon_row_size
  )

  row_sizes <- c(row_sizes, codon_row_size)

  list(
    plot = plot_result +
    ggrepel::geom_text_repel(
      data = codon_query,
      aes(
        x = genomic_position,
        label = codon
      ),
      y = 10,
      max.overlaps = Inf,
      colour = codon_query$colour,
      min.segment.length = 0,
      nudge_y = -10,
      # angle = 90,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20,
      size = 3,
      # hjust = 0,
      # segment.size = 0.2,
      # force_pull = 0, # do not pull toward data points
      # direction = "x",
      bg.color = "grey30", # shadow color
      bg.r = 0.005, # shadow radius
      arrow = grid::arrow(length = grid::unit(0.015, "npc"))
    ) +
    geom_blank(data = codon_query, aes(x = genomic_position, y = 10)) +
    ggh4x::facetted_pos_scales(
      y = list(
        grepl("annotation_plot_", name) ~ scale_y_continuous(
          limits = c(0, 10),
          breaks = NULL,
          labels = NULL
        )
      )
    ),
    row_sizes = row_sizes
  )
}

#' @importFrom ggplot2 ggplot geom_bar aes scale_color_manual coord_cartesian xlab ylab theme_minimal theme facet_grid geom_blank scale_y_continuous labeller
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggh4x facetted_pos_scales
#' @importFrom grid unit arrow
.main_plot <- function(
  track_to_plot,
  dataset_names,
  condition_names,
  plot_colours,
  region_labels,
  plot_labels,
  codon_queries,
  full_size_plots,
  half_size_plots,
  legend_position,
  text_size
)
{
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

  # print(track_to_plot)
  plot_result <- ggplot(
    track_to_plot,
    aes(y = score,
        x = genomic_position,
        colour = framing,
        fill = framing)
  ) +
    geom_bar(stat = "identity", na.rm = TRUE) +
    geom_bar(
      aes(x = genomic_position,
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
    ) +
    facet_grid(
      name ~ track_group,
      scales = "free",
      space = "free",
      labeller = labeller(track_group = region_labels, name = plot_labels),
      switch = "y"
    )

  col_size <- c(1)
  if (length(unique(track_to_plot$track_group)) == 2)
  {
    if (unique(track_to_plot$track_group)[1] == "5_UTR")
    {
      col_size <- c(0.3, 0.5)
    }
    else
    {
      col_size <- c(0.5, 0.3)
    }
  }
  else if (length(unique(track_to_plot$track_group)) == 3)
  {
    col_size <- c(0.3, 0.5, 0.3)
  }

  row_sizes <- c(rep(0.5, full_size_plots), rep(0.15, half_size_plots))
  if (!is.null(codon_queries) & nrow(codon_queries) > 0)
  {
    # split codons_queries into list based on plot_number
    codon_queries <- split(codon_queries, codon_queries$name)

    # iterate over codon queries list
    for (plot_number in seq_along(codon_queries))
    {
      # add below plot result for each element in list
      codon_query <- codon_queries[[plot_number]]
      plot_data <- .codon_plot(plot_result, codon_query, row_sizes)
      plot_result <- plot_data[[1]]

      row_sizes <- plot_data[[2]]
    }
  }

  plot_result + ggh4x::force_panelsizes(
    rows = row_sizes,
    cols = col_size
  )
}

#' @importFrom dplyr arrange mutate bind_rows distinct
.search_for_codons <- function(
  sequence,
  codons_to_search,
  annotation_label,
  plot_number,
  original_start,
  start_position,
  stop_position,
  in_frame,
  is_start
)
{
  codons_to_plot <- sapply(1:3, function(frame)
  {
    codons <- sapply(seq(from = frame, to = nchar(sequence), by = 3),
                     function(i) substr(sequence, i, i + 2))

    found <- codons %in% codons_to_search
    names(found) <- codons

    (which(found) * 3) + frame - 3
  }, simplify = F) |> unlist() |> sort()

  if (length(codons_to_plot) == 0)
    return(NULL)

  codons_to_plot <- data.frame(
    codon = names(codons_to_plot),
    genomic_position = codons_to_plot,
    p_site_framing = (codons_to_plot - 1) %% 3,
    annotation_label = annotation_label,
    plot_number = plot_number,
    name = paste0("annotation_plot_", plot_number)
  )

  codons_to_plot <- codons_to_plot |>
    dplyr::filter(genomic_position >= start_position & genomic_position <= stop_position)

  if (in_frame == T)
  {
    codons_to_plot <- codons_to_plot |>
      dplyr::filter(p_site_framing == ((original_start - 1) %% 3))
  }

  codons_to_plot |> dplyr::mutate(annotation_label = as.factor(annotation_label))
}

#' @importFrom dplyr arrange mutate bind_rows distinct
.codon_queries <- function(
  codon_queries,
  track_to_plot,
  sequence,
  start_position,
  stop_position,
  original_start,
  original_stop,
  no_orf,
  plot_colours,
  region_or_orf
)
{
  if (is.null(codon_queries) | no_orf == T)
  {
    return(data.frame())
  }

  labels <- sapply(codon_queries, `[[`, "annotation_label")
  label_to_plot <- match(labels, unique(labels))

  annotations <- lapply(seq_along(codon_queries), function(at)
  {
    query <- codon_queries[[at]]
    plot_number <- label_to_plot[at]

    codons_to_plot <- NULL

    if (is.logical(query$annotate_stop) && query$annotate_stop == T)
    {
      codons_to_plot <- .search_for_codons(
        sequence,
        c("UAA", "UAG", "UGA"),
        original_start,
        start_position,
        stop_position,
        in_frame = query$in_frame,
        is_start = F,
        annotation_label = query$annotation_label,
        plot_number = plot_number
      )
    }
    else if (is.character(query$annotate_stop))
    {
      codons_to_plot <- .search_for_codons(
        sequence,
        query$annotate_stop,
        original_start,
        start_position,
        stop_position,
        in_frame = query$in_frame,
        is_start = F,
        annotation_label = query$annotation_label,
        plot_number = plot_number
      )
    }
    else if (query$annotate_start == T)
    {
      # Create positional information for the ORFs start codon
      start_codon <- substr(sequence, start = original_start, stop = original_start + 2)
      codons_to_plot <- data.frame(
        codon = start_codon,
        genomic_position = original_start,
        plot_number = plot_number,
        is_start = T,
        annotation_label = query$annotation_label,
        name = paste0("annotation_plot_", plot_number)
      )
    }
    else if (!is.null(query$annotation_codons))
    {
      codons_to_plot <- .search_for_codons(
        sequence,
        query$annotation_codons,
        original_start,
        start_position,
        stop_position,
        in_frame = query$in_frame,
        is_start = T,
        annotation_label = query$annotation_label,
        plot_number = plot_number
      )
    }

    if (!is.null(codons_to_plot) && nrow(codons_to_plot) > 0)
    {
      if (!is.null(query$colour) && query$colour == "use_framing")
      {
        codons_to_plot$colour = plot_colours[as.character(codons_to_plot$p_site_framing)]
      }
      else if (!is.null(query$colour))
      {
        codons_to_plot$colour = query$colour
      }
      else
      {
        codons_to_plot$colour = "black"
      }
    }

    codons_to_plot
  })

  if (length(annotations) == 0)
  {
    return(data.frame())  # <- return empty, not NULL and not stop()
  }

  annotations <- Filter(Negate(is.null), annotations) |> dplyr::bind_rows()

  plot_levels <- unique(paste0("annotated_plot_", annotations$plot_number))

  annotations |>
    dplyr::arrange(genomic_position) |>
    dplyr::mutate(
      framing = "no_framing",
      name = factor(name, levels = unique(name)),
      track_group = if (is.null(original_start) & is.null(original_stop))
        "Transcript" else
          factor(
            case_when(
              !is.null(original_start) & genomic_position < original_start ~ "5_UTR",
              !is.null(original_start) & genomic_position >= original_start & genomic_position <= original_stop ~ region_or_orf,
              !is.null(original_start) & genomic_position > original_stop ~ "3_UTR"), levels = c("5_UTR", region_or_orf, "3_UTR")
          )
    ) |>
    dplyr::distinct(genomic_position, .keep_all = T)
}

.create_orf_list <- function(
    transcript_id,
    frame,
    start_codon,
    stop_codon,
    orf_start,
    orf_stop,
    nt_length
)
{
  list(
    transcript_id = transcript_id,
    frame = frame,
    start_codon = start_codon,
    stop_codon = stop_codon,
    orf_start = orf_start,
    orf_stop = orf_stop,
    nt_length = nt_length
  )
}

#' @importFrom ggplot2 ggplot geom_bar scale_color_manual aes labs theme_minimal theme scale_y_continuous facet_grid element_text element_blank labeller label_wrap_gen
#' @importFrom tidyr drop_na
.framing_plots <- function(
  track_to_plot,
  plot_colours,
  read_names,
  text_size
)
{
  track_to_plot <- track_to_plot |> drop_na()
  read_names <- read_names[names(read_names) %in% track_to_plot$name]

  track_to_plot$region <- factor(track_to_plot$region, levels = c("ORF", "Region", "Transcript", "Outside of ORF"))

  read_names[names(read_names)] <- paste(read_names[names(read_names)], "P-Sites")

  names(plot_colours) <- sub("intron_psite", "Intronic", names(plot_colours))
  track_to_plot$p_site_framing <- sub("intron_psite", "Intronic", track_to_plot$p_site_framing)

  ggplot(track_to_plot, aes(x = p_site_framing, y = p_sites, fill = p_site_framing)) +
    geom_bar(stat = "identity") +
    scale_color_manual(values = plot_colours, aesthetics = c("fill", "colour")) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      text = element_text(size = text_size),
      legend.title = element_blank(),
      legend.position = "none",
      strip.placement = "left",
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(
      name ~ region,
      # scale = "free_y",
      labeller = labeller(
        name = function(labels) label_wrap_gen(width = 15)(read_names[labels]),
        region = label_wrap_gen(width = 15)
      ),
      switch = "y"
    )
}

.integer_breaks <- function(n = 50, ...)
{
  fxn <- function(x)
  {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }

  return(fxn)
}

# internal function to remove code duplication
.check_valid_input <- function(
  orf_object,
  start_position
)
{
  if (is.null(orf_object) && is.null(start_position))
  {
    stop("Error! No ORF filter supplied, please either supply a result from inspectorORF::get_orf(), or use a transcript_filter")
  }

  orf_object_columns <- c(
    "transcript_id",
    "frame",
    "start_codon",
    "stop_codon",
    "orf_start",
    "orf_stop",
    "nt_length"
  )

  if (!is.null(orf_object) && sum(names(orf_object) %in% orf_object_columns) != 7)
  {
    stop("Error! the supplied orf details is invalid")
  }
}

.set_feature_order <- function(x)
{
  type <- ifelse(grepl("^exon", x), "exon", "intron")
  num <- as.integer(sub(".*_(\\d+)$", "\\1", x))
  order <- rank(num, ties.method = "first")
  # build combined label order interleaving exon/intron per number
  unique_nums <- sort(unique(num))
  order_labels <- as.vector(rbind(
    paste0("exon_", unique_nums),
    paste0("intron_", unique_nums)
  ))
  order_labels <- order_labels[order_labels %in% unique(x)]
  factor(x, levels = order_labels, ordered = TRUE)
}

.add_introns <- function(exon_info)
{
  exons_by_tx <- GenomicRanges::split(exon_info, exon_info$transcript_id)
  introns_by_tx <- S4Vectors::endoapply(exons_by_tx, function(exons)
  {
    reduced <- GenomicRanges::reduce(exons) |> sort()
    introns <- GenomicRanges::gaps(reduced)
    # restrict to within the transcript's bounding range
    IRanges::subsetByOverlaps(introns, range(reduced))
  })

  introns <- unlist(introns_by_tx, use.names = FALSE)
  GenomicRanges::mcols(introns)$transcript_id <- rep(names(introns_by_tx), lengths(introns_by_tx))

  introns$feature_number <- paste0("intron_", sequence(lengths(introns_by_tx)))
  # GenomicRanges::mcols(gene_info)$type <- "exon"
  # GenomicRanges::mcols(introns)$type <- "intron"

  exon_info <- c(exon_info, introns)

  exon_info <- exon_info[order(
    exon_info$transcript_id,
    GenomicRanges::seqnames(exon_info),
    ifelse(GenomicRanges::strand(exon_info) == "+", GenomicRanges::start(exon_info), -GenomicRanges::end(exon_info))
  )]

  exon_info
}

#' @importFrom ggh4x force_panelsizes
#' @importFrom dplyr filter mutate coalesce select arrange rename bind_rows case_when summarise group_by left_join
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom stats na.omit
.plot_helper <- function(
  transcript_tracks,
  orf_object,
  transcript_filter,
  start_position,
  stop_position,
  plot_region,
  plot_colours,
  scale_to_psites,
  split_exons,
  plot_transcript_summary,
  codon_queries,
  condition_names,
  plot_read_pairs,
  dataset_names,
  one_plot,
  interactive,
  legend_position,
  text_size,
  .tx_plot,
  stop_codons = c("UAG", "UAA", "UGA")
)
{
  # .check_valid_input(orf_object, start_position, plot_region)

  # if an orf has been supplied by find_orf, or get_orf, use this information instead
  # otherwise, requires transcript_filter, start position and (optionally), the stop position
  if (!is.null(orf_object))
  {
    transcript_filter <- orf_object$transcript_id
    start_position <- orf_object$orf_start
    stop_position <- orf_object$orf_stop

    if (is.null(plot_region))
    {
      plot_region = c(start_position, stop_position)
    }

    if (!(transcript_filter %in% transcript_tracks@transcript_ids))
    {
      stop(paste("Error!", transcript_filter, "not found within the BED data"))
    }
  }

  # if orf details were not supplied and only a start codon was given, find the next stop position
  if (!is.null(start_position) & is.null(stop_position))
  {
    stop_position <- .get_stop_from_start(transcript_tracks, transcript_filter, start_position, stop_codons)

    if (.tx_plot == F && is.null(plot_region[2]))
    {
      plot_region[2] = stop_position
    }
  }

  if (.tx_plot == F && (is.null(transcript_filter) || is.null(start_position) || is.null(stop_position)))
  {
    stop(paste0("Unable to obtain the transcript_id (",transcript_filter,
                 "), start position (", start_position,
                 "), or stop position (", stop_position, ")"))
  }

  original_start <- start_position
  original_stop <- stop_position

  sequence <- transcript_tracks@sequences[[transcript_filter]]
  track_to_plot <- transcript_tracks@tracks[[transcript_filter]] |> as.data.frame()

  start_position <- if (!is.na(plot_region[1]) && plot_region[1] == -1)
  {
    1
  }
  else if (!is.null(original_start) && !is.na(plot_region[1]) && original_start != plot_region[1])
  {
    original_start - plot_region[1]
  }
  else if (!is.na(plot_region[1]))
  {
    plot_region[1]
  }
  else
  {
    original_start
  }

  stop_position <- if (!is.na(plot_region[2]) && plot_region[2] == -1)
  {
    nrow(track_to_plot)
  }
  else if (!is.null(original_stop) && !is.na(plot_region[2]) && original_stop != plot_region[2])
  {
    original_stop + plot_region[2]
  }
  else if (!is.na(plot_region[2]))
  {
    plot_region[2]
  }
  else
  {
    original_stop
  }

  no_orf <- F
  if (is.null(original_start) & is.null(original_stop))
  {
    plot_transcript_summary <- F
    no_orf <- T
    original_start <- start_position
    original_stop <- stop_position
  }

  found_read_names <- track_to_plot$name |> unique()
  read_names_count <- found_read_names |> length()

  if (.tx_plot & no_orf)
  {
    orf_or_region <- "Region"

    if (is.character(split_exons) & split_exons == "with_introns")
    {
      print("Adding introns to plot")
      track_to_plot <- track_to_plot |> dplyr::mutate(track_group = feature_number)

      # Extract unique region names
      groups <- unique(track_to_plot$track_group)

      # Assign labels: "exon" if name starts with exon_, else "intron"
      region_labels <- setNames(
        ifelse(grepl("^exon", groups), "Exon", "Intron"),
        groups
      )
    }
    else if (is.logical(split_exons) & split_exons == T)
    {
      track_to_plot <- track_to_plot |> dplyr::filter(grepl("^exon", feature_number)) |>
        dplyr::mutate(track_group = feature_number)

      region_labels <- setNames(
        rep("Exon", length(unique(track_to_plot$track_group))),
        unique(track_to_plot$track_group)
      )
    }
    else
    {
      region_labels <- c("5_UTR" = "5'", "Region" = "Region", "3_UTR" = "3'")

      track_to_plot <- track_to_plot |> dplyr::filter(grepl("^exon", feature_number)) |>
        dplyr::mutate(track_group = "Transcript")
    }
  }
  else
  {
    region_labels <- c("5_UTR" = "5'", "ORF" = "ORF", "3_UTR" = "3'")
    orf_or_region <- "ORF"

    track_to_plot <- track_to_plot |> dplyr::filter(grepl("^exon", feature_number)) |>
      dplyr::mutate(
        track_group = factor(
          dplyr::case_when(
            genomic_position < original_start ~ "5_UTR",
            genomic_position >= original_start & genomic_position <= original_stop ~ orf_or_region,
            genomic_position > original_stop ~ "3_UTR"), levels = c("5_UTR", orf_or_region, "3_UTR")
        )
      )
  }

  group_levels <- track_to_plot |>
    dplyr::distinct(track_group, genomic_position) |>
    dplyr::arrange(genomic_position) |>
    dplyr::pull(track_group) |>
    unique()

  track_to_plot <- track_to_plot |>
    dplyr::mutate(
      track_group = factor(track_group, levels = group_levels)
    )

  filtered_framed_tracks <- intersect(transcript_tracks@framed_tracks, names(plot_read_pairs))

  if (!is.null(dataset_names))
  {
    track_to_plot <- track_to_plot |> dplyr::filter(name %in% c(filtered_framed_tracks, names(dataset_names)))
  }

  # track_to_plot <- track_to_plot |> arrange(seqnames, start, strand, genomic_position, framing)

  p_sites <- track_to_plot |> dplyr::select(transcript_id, seqnames, start, end, name, score, framing) |>
    dplyr::filter(name %in% names(plot_read_pairs)) |>
    dplyr::rename(p_sites = score,
                  p_site_framing = framing) |>
    dplyr::filter(!is.na(p_site_framing)) |>
    dplyr::arrange(match(name, names(plot_read_pairs))) |>
    dplyr::mutate(name = rep(plot_read_pairs, each = length(name) / length(plot_read_pairs)))

  # print(track_to_plot)
  track_to_plot <- track_to_plot |> dplyr::left_join(p_sites, by = c("transcript_id", "seqnames", "start", "end", "name")) |>
    dplyr::filter(!(name %in% names(plot_read_pairs))) |>
    dplyr::mutate(p_site_framing = dplyr::coalesce(p_site_framing, name))

  if (any(!unique(track_to_plot$name) %in% names(plot_colours)) == T)
  {
    names_to_change <- unique(track_to_plot$name)[!(unique(track_to_plot$name) %in% names(plot_colours))]
    track_to_plot[track_to_plot$name == names_to_change, "framing"] <- track_to_plot[track_to_plot$name == names_to_change, "og_framing"]
  }

  # print(track_to_plot)

  plot_labels <- dataset_names[!(names(dataset_names) %in% filtered_framed_tracks)]
  plot_labels[names(plot_labels)] <- ifelse(
    names(plot_labels) %in% names(condition_names) & condition_names[names(plot_labels)] != "",
    paste0(condition_names[names(plot_labels)], "\n", dataset_names[names(plot_labels)]),
    paste0(dataset_names[names(plot_labels)])
  )

  plot_labels[plot_read_pairs] <- paste0(plot_labels[plot_read_pairs], "\n+", dataset_names[names(plot_read_pairs)])

  track_to_plot <- track_to_plot |> dplyr::mutate(name = factor(name, levels = names(plot_labels)))

  codon_queries <- .codon_queries(
    codon_queries,
    track_to_plot,
    sequence,
    start_position,
    stop_position,
    original_start,
    original_stop,
    no_orf,
    plot_colours,
    orf_or_region
  )

  if (nrow(codon_queries) > 0)
  {
    codon_labels <- setNames(unique(codon_queries$annotation_label), unique(codon_queries$name))
    plot_labels <- c(plot_labels, codon_labels)
  }

  # obtain the tracks to plot
  to_plot <- track_to_plot |> dplyr::filter(genomic_position >= start_position & genomic_position <= stop_position)

  if (scale_to_psites == T)
  {
    y_limit <- to_plot |> dplyr::filter(name %in% plot_read_pairs) |>
      dplyr::group_by(name) |>
      dplyr::summarise(max_p = max(p_sites |> na.omit()))

    to_plot <- to_plot |>
      dplyr::left_join(y_limit, by = "name") |>
      dplyr::mutate(max_p = coalesce(max_p, score)) |>
      (\(df) dplyr::mutate(df, score = pmin(df$max_p, df$score)))()
  }

  full_size_plots <- length(plot_read_pairs)
  paired_datasets <- union(names(plot_read_pairs), unname(plot_read_pairs))
  half_size_plots <- length(setdiff(names(dataset_names), paired_datasets))

  plot_result <- .main_plot(
    to_plot,
    dataset_names,
    condition_names,
    plot_colours,
    region_labels,
    plot_labels,
    codon_queries,
    full_size_plots,
    half_size_plots,
    legend_position,
    text_size
  ) + theme(legend.position = legend_position)

  orf_tracks <- track_to_plot |>
    dplyr::filter(genomic_position >= original_start & genomic_position <= original_stop) |>
    dplyr::mutate(region = ifelse(.tx_plot == TRUE & no_orf, "Transcript", orf_or_region))

  if (plot_transcript_summary)
  {
    rest_of_transcript <- track_to_plot |>
      dplyr::filter(!(genomic_position %in% orf_tracks$genomic_position)) |>
      dplyr::mutate(region = "Outside of ORF")

    combined_tracks <- dplyr::bind_rows(orf_tracks, rest_of_transcript)
  }
  else
  {
    combined_tracks <- orf_tracks
  }

  orf_frame_plot <- .framing_plots(
    combined_tracks,
    plot_colours,
    condition_names,
    text_size
  ) + theme(legend.position = "none")

  plot_list <- list(plot_result, orf_frame_plot)

  plot_figure_width <- c(3, 0.75)

  if (interactive && requireNamespace("plotly", quietly = TRUE))
  {
    plot_list <- lapply(plot_list, plotly::ggplotly)
  }
  else if (interactive)
  {
    message("Plotly not installed - falling back to static figure.")
  }

  if (one_plot)
  {
    if (interactive && requireNamespace("plotly", quietly = TRUE))
    {
      return(
        plotly::subplot(
          plot_list,
          nrows = 1,
          ncols = length(plot_list),
          shareX = FALSE,
          shareY = FALSE,
          titleX = TRUE,
          titleY = TRUE
        )
      )
    }
    else
    {
      patchwork_plot <- patchwork::wrap_plots(plot_list, nrow = 1, widths = plot_figure_width) +
        # patchwork::plot_layout(guides = "collect") #+
        patchwork::plot_annotation()
        # theme(legend.position = "bottom")

      return(patchwork_plot)
    }
  }

  plot_list
}
