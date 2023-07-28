
pacman::p_load(magrittr, tidyverse, UWCCC.GWAS)

pacman::p_load(
  snpStats, ggbio,
  EnsDb.Hsapiens.v75, AnnotationDbi,
  org.Hs.eg.db,
  TxDb.Hsapiens.UCSC.hg19.knownGene,
  BSgenome.Hsapiens.UCSC.hg19,
  biomaRt)


# get significant snps
vitd_assoc <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "vitd_associations.qs"))


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- org.Hs.eg.db 
ensdb <- EnsDb.Hsapiens.v75

snp_dist <- 250e3
pval_thr <- 1e-8

future::plan(future::multiprocess, workers = 20)

sccs_vdbp_loci <-
  qs::qread(here::here("results", "2020_10_20_vitd_fine_mapping", "qs",
    str_c("sccs_vitdbp_loci_snpdist",
      scales::comma(snp_dist, scale = 1e-3, suffix = "kp"), "_pval_thr",
      scales::scientific(pval_thr),
    ".qs")))

gwas_dir <- here::here("extracted_data", "SCCS", "imputation", 
  "merge_impute", "qc_filtered")

files <- list.files(gwas_dir, pattern = "rsq50", full.names = TRUE)

get_plink_files <- function(all_files, chrom) {

  dir <- dirname(all_files) %>%
    unique()
  all_files <- basename(all_files) %>%
    str_subset(str_c("chr", chrom, "_"))
  file.path(dir, all_files)
}

# this gets genotype for every loci
sccs_vdbp_loci %<>%
  mutate(
    plink_files = map(chrom, ~ get_plink_files(files, .)),
    snpmat = map2(plink_files, locus,
     ~ snpStats::read.plink(
        bed = str_subset(.x, ".bed"),
        bim = str_subset(.x, ".bim"),
        fam = str_subset(.x, ".fam"),
        select.snps = .y$id)))

compute_ld <- function(locus, snpmat) {

  ld <- cor(as(snpmat$genotype, "numeric"))
  index_snp <- locus %>%
    dplyr::filter(index) %>%
    dplyr::pull(id)

  ld <- ld[, index_snp]
  locus %>%
    dplyr::mutate(ld = ld[locus$id])

}


get_locus_range <- function(chrom, locus) {

  GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges(
      start = min(locus$pos),
      end = max(locus$pos)))

}


locus_zoom_plot <- function(chrom, locus, db) {

  range <- get_locus_range(chrom, locus)
  gene_annot <- ggbio::autoplot(db,
      AnnotationFilter::GRangesFilter(range),
      names.expr = "gene_name")
  main_plot <- locus %>%
    ggplot(aes(x = pos, y = -log10(p), fill = ld^2)) +
      geom_point(aes(shape = index), size = 5, colour = "black") +
      scale_shape_manual(values = c(21, 23), guide = FALSE) +
      labs(fill = expression(r^2), y = expression(-log[10](p.value))) +
      geom_hline(yintercept = 8, linetype = 2) +
      scale_fill_distiller(palette = "RdYlGn",
        labels = scales::percent_format(1))

  tracks(main_plot, gene_annot, heights = c(3, 1), xlim = range,
    theme = theme_bw() + theme(
        legend.key.width = unit(.05, "npc"),
        legend.position = "top")) +
    scale_x_sequnit("kb")
}

pngdir <- here::here("results", "zz_methods", "figs")

sccs_vdbp_loci %<>%
  dplyr::mutate(
    locus = map(locus, ~ mutate(., index = pos %in% vitd_assoc[[1]]$pos))) %>%
  dplyr::filter(map_lgl(locus, ~ any(.$index))) %>%
  dplyr::mutate(
    locus = furrr::future_map2(locus, snpmat, compute_ld),
    locus_plot = map2(chrom, locus, locus_zoom_plot,
      ensdb),
    min = map_int(locus, ~ min(.$pos)),
    max = map_int(locus, ~ max(.$pos)),
    min = scales::comma(min, accuracy = 1, scale = 1e-3, suffix = "kb"),
    max = scales::comma(max, accuracy = 1, scale = 1e-3, suffix = "kb"),
    filename = glue::glue("{dr}/SCCS_vdbp_{chrom}_{min}_{max}_locuszoom.png",
      dr = pngdir),
    walk2(filename, locus_plot, ggsave,
      width = 10, height = 6, units = "in"))
