#!/usr/bin/env Rscript

###############################################################################
# run_qdnaseq_genome_integrity.R
#
# Purpose:
#   Run a QDNAseq-based genome integrity screen from BAM files and produce
#   genome-wide copy-number plots suitable for cultured-cell QC
#   (e.g. iPSCs, fibroblasts, ESCs).
#
# Main outputs per sample:
#   - genomewide_<sample>.jpg
#   - <sample>.qdnaseq.rds
#   - <sample>_bins.tsv.gz
#   - <sample>_segments.tsv
#
# Run-level outputs:
#   - sample_qc_summary.tsv
#   - SESSION_INFO.txt
#
# Notes:
#   - Copy number is plotted as log2 ratio:
#       ~0   = copy-number neutral / diploid baseline
#       >0   = relative gain
#       <0   = relative loss
#   - Intended for broad genome integrity QC, not focal CNV discovery.
#
# Run from RStudio's Terminal tab (or any shell) with:
#   Rscript run_qdnaseq_genome_integrity.R 
#     --bam-dir bam
#     --out-dir qdnaseq_output
#     --bin-size 100
#     --genome hg38
#     --sex-map sample_sex.csv
#     --map-cutoff 80
#     --bases-cutoff 99.9
#     --seed 170826
#
###############################################################################

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(optparse)
  library(tools)
  library(grid)
})

SCRIPT_VERSION <- "1.0.1"

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------

option_list <- list(
  make_option(
    "--bam-dir",
    dest = "bam_dir",
    type = "character",
    default = "bam",
    help = "Directory containing BAM (.bam) and index (.bai) files [default: %default]"
  ),
  make_option(
    "--out-dir",
    dest = "out_dir",
    type = "character",
    default = "qdnaseq_output",
    help = "Output directory [default: %default]"
  ),
  make_option(
    "--bin-size",
    dest = "bin_size",
    type = "integer",
    default = 100,
    help = "Bin size in kb for QDNAseq [default: %default]"
  ),
  make_option(
    "--genome",
    dest = "genome",
    type = "character",
    default = "hg38",
    help = "Genome build used for bin annotations, e.g. hg19 or hg38 [default: %default]"
  ),
  make_option(
    "--sex-map",
    dest = "sex_map",
    type = "character",
    default = NULL,
    help = "REQUIRED. CSV with columns sample,sex where sex is XX or XY. Every sample in --bam-dir must have an entry."
  ),
  make_option(
    "--include-mt",
    dest = "include_mt",
    action = "store_true",
    default = FALSE,
    help = "Include mitochondrial chromosome MT [default: %default]"
  ),
  make_option(
    "--no-sex-ploidy-correction",
    dest = "no_sex_ploidy_correction",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Disable the XY sex-chromosome ploidy correction (rescales chrX/chrY",
      "copy number by 2x for XY samples so they sit at the same baseline",
      "as autosomes). The correction is applied by default; set this flag",
      "to see the raw, uncorrected QDNAseq output instead [default: %default]"
    )
  ),
  make_option(
    "--map-cutoff",
    dest = "map_cutoff",
    type = "double",
    default = 80,
    help = "Minimum mappability percentage for QC flagging [default: %default]"
  ),
  make_option(
    "--bases-cutoff",
    dest = "bases_cutoff",
    type = "double",
    default = 99.9,
    help = "Minimum characterized bases percentage for QC flagging [default: %default]"
  ),
  make_option(
    "--min-mapq",
    dest = "min_mapq",
    type = "integer",
    default = 37,
    help = "Minimum BAM mapping quality passed to binReadCounts() [default: %default]"
  ),
  make_option(
    "--segment-alpha",
    dest = "segment_alpha",
    type = "double",
    default = 0.001,
    help = "DNAcopy alpha used by segmentBins() [default: %default]"
  ),
  make_option(
    "--segment-undo-sd",
    dest = "segment_undo_sd",
    type = "double",
    default = 0.5,
    help = "DNAcopy undo.SD used by segmentBins() [default: %default]"
  ),
  make_option(
    "--segment-tol",
    dest = "segment_tol",
    type = "double",
    default = 0.01,
    help = "Tolerance for merging adjacent equal segment values in plotting [default: %default]"
  ),
  make_option(
    "--seed",
    dest = "seed",
    type = "integer",
    default = 170826,
    help = "Random seed for reproducible DNAcopy segmentation [default: %default]"
  ),
  make_option(
    "--plot-filtered-bins",
    dest = "plot_filtered_bins",
    type = "logical",
    default = TRUE,
    help = "Show filtered (low-quality) bins as light grey points: TRUE or FALSE [default: %default]"
  ),
  make_option(
    "--width-mm",
    dest = "width_mm",
    type = "double",
    default = 95,
    help = "Plot width in mm [default: %default]"
  ),
  make_option(
    "--height-mm",
    dest = "height_mm",
    type = "double",
    default = 60,
    help = "Plot height in mm [default: %default]"
  ),
  make_option(
    "--ymin",
    dest = "ymin",
    type = "double",
    default = -1.2,
    help = "Plot lower y-limit [default: %default]"
  ),
  make_option(
    "--ymax",
    dest = "ymax",
    type = "double",
    default = 1.2,
    help = "Plot upper y-limit [default: %default]"
  ),
  make_option(
    "--verbose",
    dest = "verbose",
    type = "logical",
    default = TRUE,
    help = "Print progress messages: TRUE or FALSE [default: %default]"
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list
)

opt <- parse_args(parser)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

msg <- function(...) {
  if (isTRUE(opt$verbose)) {
    message(...)
  }
}

stop_if_not_dir <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !dir.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

validate_writable_outdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  testfile <- file.path(path, paste0(".write_test_", Sys.getpid()))
  ok <- tryCatch(
    {
      file.create(testfile)
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  if (!isTRUE(ok)) {
    stop("Output directory is not writable: ", path, call. = FALSE)
  }
  unlink(testfile)
}

read_sex_map <- function(path) {
  if (!file.exists(path)) {
    stop("Sex map file not found: ", path, call. = FALSE)
  }
  
  x <- read.csv(path, stringsAsFactors = FALSE)
  
  if (!all(c("sample", "sex") %in% names(x))) {
    stop("Sex map must contain columns: sample, sex", call. = FALSE)
  }
  
  x$sample <- trimws(gsub("^\ufeff", "", as.character(x$sample)))
  x$sex <- toupper(trimws(as.character(x$sex)))
  
  bad <- setdiff(unique(x$sex), c("XX", "XY"))
  if (length(bad) > 0) {
    stop("Invalid sex values in sex map: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  
  if (anyDuplicated(x$sample)) {
    dup <- unique(x$sample[duplicated(x$sample)])
    stop("Duplicate sample names in sex map: ", paste(dup, collapse = ", "), call. = FALSE)
  }
  
  x
}

safe_log2 <- function(x, min_positive = 1e-3) {
  ifelse(is.finite(x) & x >= min_positive, log2(x), NA_real_)
}

chr_levels <- function(chr_vec, include_Y = TRUE, include_MT = FALSE) {
  autos <- as.character(1:22)
  lev <- c(
    autos,
    "X",
    if (isTRUE(include_Y)) "Y" else NULL,
    if (isTRUE(include_MT)) "MT" else NULL
  )
  factor(as.character(chr_vec), levels = lev[lev %in% unique(as.character(chr_vec))])
}

mat_to_long <- function(mat, value_name) {
  df <- cbind(bin = seq_len(nrow(mat)), as.data.frame(mat, check.names = FALSE))
  pivot_longer(df, -bin, names_to = "sample", values_to = value_name)
}

apply_sex_ploidy_correction <- function(df, is_xy = FALSE) {
  if (!isTRUE(is_xy)) {
    return(df)
  }
  df |>
    mutate(
      CN_raw = ifelse(chromosome %in% c("X", "Y"), CN_raw * 2, CN_raw),
      SEG_raw = ifelse(chromosome %in% c("X", "Y"), SEG_raw * 2, SEG_raw)
    )
}

extract_plot_df <- function(cn, include_Y = TRUE, include_MT = FALSE, is_xy = include_Y) {
  fd <- as.data.frame(Biobase::fData(cn))
  fd$bin <- seq_len(nrow(fd))
  
  bases_pct <- as.numeric(fd$bases)
  map_raw <- as.numeric(fd$mappability)
  
  mx <- suppressWarnings(max(map_raw, na.rm = TRUE))
  mappability_pct <- if (is.finite(mx) && mx <= 1.1) 100 * map_raw else map_raw
  mappability_pct <- pmin(pmax(mappability_pct, 0), 100)
  
  bins <- tibble(
    bin = fd$bin,
    chromosome = chr_levels(fd$chromosome, include_Y = include_Y, include_MT = include_MT),
    start = as.numeric(fd$start),
    end = as.numeric(fd$end),
    bases_pct = bases_pct,
    mappability_pct = mappability_pct
  ) |>
    filter(!is.na(chromosome))
  
  assays <- Biobase::assayDataElementNames(cn)
  
  cn_long <- mat_to_long(Biobase::assayDataElement(cn, "copynumber"), "CN_raw")
  
  seg_long <- if ("segmented" %in% assays) {
    mat_to_long(Biobase::assayDataElement(cn, "segmented"), "SEG_raw")
  } else {
    tibble(bin = integer(), sample = character(), SEG_raw = numeric())
  }
  
  out <- bins |>
    left_join(cn_long, by = "bin") |>
    left_join(seg_long, by = c("bin", "sample")) |>
    arrange(sample, chromosome, start) |>
    as_tibble() |>
    apply_sex_ploidy_correction(is_xy = is_xy)

  out |>
    mutate(
      CN_log2  = safe_log2(CN_raw),
      SEG_log2 = safe_log2(SEG_raw)
    )
}

add_qc_flags <- function(df, bases_cut, map_cut) {
  df |>
    mutate(
      low_bases = is.finite(bases_pct) & bases_pct < bases_cut,
      low_map = is.finite(mappability_pct) & mappability_pct < map_cut,
      status = case_when(
        low_bases & low_map ~ "Both low",
        low_bases ~ "Bases < cutoff",
        low_map ~ "Low mappability",
        TRUE ~ "OK"
      )
    )
}

get_chr_offsets <- function(cn, include_Y = TRUE, include_MT = FALSE) {
  fa <- as.data.frame(Biobase::fData(cn))
  fa$chromosome <- chr_levels(fa$chromosome, include_Y = include_Y, include_MT = include_MT)
  fa <- subset(fa, !is.na(chromosome))
  fa$end <- as.numeric(fa$end)
  
  fa |>
    group_by(chromosome) |>
    summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") |>
    arrange(chromosome) |>
    mutate(
      chr_len = as.numeric(chr_len),
      offset = lag(cumsum(chr_len), default = 0),
      xmin_mb = offset / 1e6,
      xmax_mb = (offset + chr_len) / 1e6,
      center_mb = (offset + chr_len / 2) / 1e6,
      band = (row_number() %% 2) == 1
    )
}

add_pos_cum <- function(df_in, offs) {
  df_in |>
    mutate(
      chromosome = factor(as.character(chromosome), levels = levels(offs$chromosome)),
      start = as.numeric(start),
      end = as.numeric(end)
    ) |>
    left_join(select(offs, chromosome, offset), by = "chromosome") |>
    mutate(
      pos_cum = (start + offset) / 1e6,
      pos_cum_end = (end + offset) / 1e6
    ) |>
    arrange(chromosome, start)
}

compute_segment_runs <- function(seg_df, tol = 1e-3) {
  seg_df |>
    filter(
      is.finite(SEG_log2),
      is.finite(pos_cum),
      is.finite(pos_cum_end)
    ) |>
    arrange(chromosome, start) |>
    group_by(chromosome) |>
    mutate(
      new_run = row_number() == 1 | abs(SEG_log2 - lag(SEG_log2)) > tol
    ) |>
    mutate(run_id = cumsum(new_run)) |>
    group_by(chromosome, run_id) |>
    summarise(
      seg_y = first(SEG_log2),
      x_minMb = first(pos_cum),
      x_maxMb = last(pos_cum_end),
      .groups = "drop"
    )
}

limited_chr_labels <- function(offs) {
  wanted <- c(as.character(1:10), "12", "14", "16", "18", "21", "X", "Y")
  labs <- as.character(offs$chromosome)
  labs[!(labs %in% wanted)] <- ""
  list(breaks = offs$center_mb, labels = labs)
}

plot_genomewide <- function(ok_df, filtered_df, seg_runs, offs, ymin = -1.2, ymax = 1.2,
                            plot_filtered_bins = FALSE) {
  total_len_mb <- max(offs$xmax_mb, na.rm = TRUE)
  lab <- limited_chr_labels(offs)
  chr_boundaries <- offs$xmax_mb[-nrow(offs)]
  
  p <- ggplot() +
    geom_rect(
      data = offs,
      aes(xmin = xmin_mb, xmax = xmax_mb, ymin = -Inf, ymax = Inf, fill = band),
      inherit.aes = FALSE,
      alpha = 0.05,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = chr_boundaries,
      linewidth = 0.25,
      colour = "grey70"
    ) +
    geom_point(
      data = ok_df,
      aes(x = pos_cum, y = CN_log2),
      colour = "black",
      alpha = 1,
      size = 0.40,
      shape = 16,
      stroke = 0
    )
  
  if (isTRUE(plot_filtered_bins) && nrow(filtered_df) > 0) {
    p <- p +
      geom_point(
        data = filtered_df,
        aes(x = pos_cum, y = CN_log2),
        colour = "grey70",
        alpha = 0.4,
        size = 0.35,
        shape = 16,
        stroke = 0
      )
  }
  
  p +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "seagreen3"
    ) +
    geom_segment(
      data = seg_runs,
      aes(x = x_minMb, xend = x_maxMb, y = seg_y, yend = seg_y),
      colour = "#6A0DAD",
      linewidth = 1.6,
      alpha = 0.98,
      lineend = "butt"
    ) +
    scale_fill_manual(values = c("TRUE" = "grey80", "FALSE" = "white")) +
    scale_y_continuous(
      "Relative copy number (log2 ratio)",
      breaks = seq(-2, 2, by = 0.5),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_x_continuous(
      "Chromosome",
      breaks = lab$breaks,
      labels = lab$labels,
      expand = expansion(mult = c(0, 0)),
      limits = c(0, total_len_mb)
    ) +
    coord_cartesian(ylim = c(ymin, ymax), clip = "on") +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks.length = unit(2, "pt"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 9, margin = margin(t = 4)),
      axis.title.y = element_text(size = 9, margin = margin(r = 6)),
      plot.margin = margin(3, 6, 2, 6)
    )
}

sample_qc_summary <- function(df) {
  tibble(
    n_bins = nrow(df),
    n_finite = sum(is.finite(df$CN_log2)),
    frac_negative = mean(df$CN_log2 < 0, na.rm = TRUE),
    median_log2 = median(df$CN_log2, na.rm = TRUE),
    median_bases_pct = median(df$bases_pct, na.rm = TRUE),
    median_mappability_pct = median(df$mappability_pct, na.rm = TRUE),
    frac_low_bases = mean(df$low_bases, na.rm = TRUE),
    frac_low_map = mean(df$low_map, na.rm = TRUE)
  )
}

write_session_info <- function(out_dir, seed) {
  lines <- c(
    paste("Script version:", SCRIPT_VERSION),
    paste("Date:", as.character(Sys.time())),
    paste("Command:", paste(commandArgs(), collapse = " ")),
    paste("Random seed:", seed),
    "",
    capture.output(sessionInfo())
  )
  writeLines(lines, con = file.path(out_dir, "SESSION_INFO.txt"))
}

validate_inputs <- function(bamfiles, bamnames, sex_map) {
  if (!length(bamfiles)) {
    stop("No BAM files found.", call. = FALSE)
  }
  
  if (anyDuplicated(bamnames)) {
    dup <- unique(bamnames[duplicated(bamnames)])
    stop(
      "Duplicate sample names inferred from BAM filenames: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }
  
  missing_bai <- bamfiles[!file.exists(paste0(bamfiles, ".bai"))]
  if (length(missing_bai)) {
    stop(
      "Missing BAM index (.bai) for: ",
      paste(basename(missing_bai), collapse = ", "),
      call. = FALSE
    )
  }
  
  not_found <- setdiff(sex_map$sample, bamnames)
  if (length(not_found)) {
    stop(
      "Sample(s) in sex map do not match any BAM-derived sample name: ",
      paste(not_found, collapse = ", "), ". ",
      "BAM-derived sample names are: ", paste(bamnames, collapse = ", "), ". ",
      call. = FALSE
    )
  }
  
  missing_from_map <- setdiff(bamnames, sex_map$sample)
  if (length(missing_from_map)) {
    stop(
      "Sample(s) found in BAM directory but missing from sex map: ",
      paste(missing_from_map, collapse = ", "), ". ",
      "Every BAM sample must have a sex assignment in the sex-map CSV.",
      call. = FALSE
    )
  }
}

export_table_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---------------------------------------------------------------------------
# Input discovery and validation
# ---------------------------------------------------------------------------

stop_if_not_dir(opt$bam_dir, "BAM directory")
validate_writable_outdir(opt$out_dir)

if (is.null(opt$sex_map) || !nzchar(opt$sex_map)) {
  stop(
    "--sex-map is required. Provide a CSV with columns sample,sex (XX or XY) ",
    "covering every sample in --bam-dir.",
    call. = FALSE
  )
}

sex_map <- read_sex_map(opt$sex_map)

bam_dir_abs <- normalizePath(opt$bam_dir, winslash = "/", mustWork = TRUE)
out_dir_abs <- normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE)

bamfiles <- list.files(
  path = bam_dir_abs,
  pattern = "\\.bam$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)
bamfiles <- sort(bamfiles)

bamnames <- file_path_sans_ext(basename(bamfiles))

validate_inputs(bamfiles, bamnames, sex_map)

msg("Found ", length(bamfiles), " BAM file(s).")
msg("BAM directory: ", bam_dir_abs)
msg("Samples: ", paste(bamnames, collapse = ", "))
msg("Genome: ", opt$genome, " | Bin size: ", opt$bin_size, " kb")
msg("Output directory: ", out_dir_abs)

# ---------------------------------------------------------------------------
# Load bin annotations
# ---------------------------------------------------------------------------

msg("Loading bin annotations: genome=", opt$genome, ", bin-size=", opt$bin_size, " kb")
bins <- tryCatch(
  getBinAnnotations(binSize = opt$bin_size, genome = opt$genome),
  error = function(e) {
    stop(
      "Could not load QDNAseq bin annotations for genome='", opt$genome,
      "', bin-size=", opt$bin_size, " kb. ",
      "Make sure the matching QDNAseq annotation package is installed.",
      call. = FALSE
    )
  }
)

# ---------------------------------------------------------------------------
# Per-sample QDNAseq processing, plotting and export
# ---------------------------------------------------------------------------

set.seed(opt$seed)
msg("Random seed set to ", opt$seed, " (DNAcopy CBS uses permutation testing).")

qc_rows <- list()

for (i in seq_along(bamfiles)) {
  bamfile <- bamfiles[i]
  sample_name <- bamnames[i]
  
  msg("Processing sample: ", sample_name)
  msg("  BAM: ", bamfile)
  
  include_Y <- sex_map$sex[match(sample_name, sex_map$sample)] == "XY"
  
  chrom_filter <- if (include_Y) NA else "Y"
  
  rc <- binReadCounts(
    bins = bins,
    bamfiles = bamfile,
    bamnames = sample_name,
    minMapq = opt$min_mapq,
    verbose = opt$verbose
  )
  
  rc <- applyFilters(
    rc,
    residual = TRUE,
    blacklist = TRUE,
    mappability = opt$map_cutoff,
    bases = opt$bases_cutoff,
    chromosomes = chrom_filter,
    verbose = opt$verbose
  )
  
  cn <- estimateCorrection(rc)
  cn <- correctBins(cn)
  cn <- normalizeBins(cn)
  cn <- smoothOutlierBins(cn)
  cn <- segmentBins(
    cn,
    alpha = opt$segment_alpha,
    undo.splits = "sdundo",
    undo.SD = opt$segment_undo_sd,
    verbose = opt$verbose
  )
  
  apply_ploidy_correction <- include_Y && !isTRUE(opt$no_sex_ploidy_correction)
  
  plot_df <- extract_plot_df(
    cn,
    include_Y = include_Y,
    include_MT = isTRUE(opt$include_mt),
    is_xy = apply_ploidy_correction
  ) |>
    filter(sample == sample_name) |>
    add_qc_flags(
      bases_cut = opt$bases_cutoff,
      map_cut = opt$map_cutoff
    )
  
  offs <- get_chr_offsets(
    cn,
    include_Y = include_Y,
    include_MT = isTRUE(opt$include_mt)
  )
  
  plot_df_pos <- add_pos_cum(plot_df, offs)
  
  ok_bins <- plot_df_pos |>
    filter(
      status == "OK",
      is.finite(CN_log2),
      is.finite(pos_cum)
    )
  
  filtered_bins <- plot_df_pos |>
    filter(
      status != "OK",
      is.finite(CN_log2),
      is.finite(pos_cum)
    )
  
  seg_input <- plot_df_pos |>
    filter(
      is.finite(SEG_log2),
      is.finite(pos_cum),
      is.finite(pos_cum_end)
    )
  
  seg_runs <- compute_segment_runs(seg_input, tol = opt$segment_tol)
  
  p <- plot_genomewide(
    ok_df = ok_bins,
    filtered_df = filtered_bins,
    seg_runs = seg_runs,
    offs = offs,
    ymin = opt$ymin,
    ymax = opt$ymax,
    plot_filtered_bins = isTRUE(opt$plot_filtered_bins)
  )
  
  out_plot <- file.path(opt$out_dir, paste0("genomewide_", make.names(sample_name), ".jpg"))
  grDevices::jpeg(
    filename = out_plot,
    width = opt$width_mm / 25.4,
    height = opt$height_mm / 25.4,
    units = "in",
    res = 600,
    quality = 100,
    bg = "white"
  )
  print(p)
  grDevices::dev.off()
  
  out_rds <- file.path(opt$out_dir, paste0(make.names(sample_name), ".qdnaseq.rds"))
  saveRDS(cn, out_rds)
  
  out_bins <- file.path(opt$out_dir, paste0(make.names(sample_name), "_bins.tsv.gz"))
  export_table_gz(
    plot_df_pos |>
      select(
        sample, chromosome, start, end, bases_pct, mappability_pct,
        low_bases, low_map, status, CN_raw, SEG_raw, CN_log2, SEG_log2,
        pos_cum, pos_cum_end
      ),
    out_bins
  )
  
  out_segments <- file.path(opt$out_dir, paste0(make.names(sample_name), "_segments.tsv"))
  write.table(seg_runs, file = out_segments, sep = "\t", row.names = FALSE, quote = FALSE)
  
  qc_rows[[sample_name]] <- sample_qc_summary(plot_df_pos) |>
    mutate(
      sample = sample_name,
      sex_mode = if (include_Y) "XY" else "XX",
      sex_ploidy_corrected = apply_ploidy_correction,
      plot_file = basename(out_plot),
      qdnaseq_rds = basename(out_rds),
      bins_tsv_gz = basename(out_bins),
      segments_tsv = basename(out_segments)
    )
  
  msg("  Wrote plot: ", out_plot)
  msg("  Wrote object: ", out_rds)
  msg("  Wrote bins: ", out_bins)
  msg("  Wrote segments: ", out_segments)
}

qc_summary <- bind_rows(qc_rows) |>
  select(sample, everything())

write.table(
  qc_summary,
  file = file.path(opt$out_dir, "sample_qc_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write_session_info(opt$out_dir, opt$seed)

msg("Done.")

