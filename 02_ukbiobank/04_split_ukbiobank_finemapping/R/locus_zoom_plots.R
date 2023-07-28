
# rs10805046	NPFFR2	4	72,888,168	G/A	0.420	-0.183	0.029	3.78E-10	1.53%
# rs10938008	NPFFR2	4	73,065,084	A/C	0.224	0.208	0.034	7.60E-10	1.47%
library(magrittr)
library(tidyverse)

source(here::here("src/gwas_funs_redo.R"))

# load summary statistics
stats0 <- ukb_afroc_gwas <-
  here::here("results/2022_08_04_ukbiobank_finemapping/qs",
   "afroc_ukb_stats_filter.qs") %>%
  qs::qread()

# fine-mapping parameters
pval <- 5e-8
window <- 250e3

stats1 <- stats0 %>%
  dplyr::filter(p <= pval)

# > stats1 %>% dplyr::pull(pos) %>% range
# [1] 72608364 72681824
# > stats1 %>% dplyr::pull(pos) %>% range %>% diff
# [1] 73460
# > stats1 %>% dplyr::pull(pos) %>% summary
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 72608364 72617612 72618049 72631222 72637197 72681824 

top_snp <- stats1 %>%
  dplyr::group_by(chrom) %>%
  dplyr::top_n(1, wt = -p) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    chrom = str_c("0", chrom))


# top_snp %<>%
#   dplyr::bind_rows(
#     top_snp %>%
#       dplyr::select(chrom, pos) %>%
#       dplyr::inner_join(
#         stats1, by = "chrom", suffix = c("ref", "")) %>%
#       dplyr::filter(abs(posref - pos) > window / 2) %>%
#       dplyr::select(-posref) %>%
#       dplyr::group_by(chrom) %>%
#       dplyr::top_n(1, wt = -p))

# Bioconductor part
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)
library(Matrix)

top_snp %<>%
  dplyr::mutate(id2 = id) %>%
  tidyr::nest(data = -c(id2)) %>%
  dplyr::mutate(
    gr = purrr::map(data, create_gr),
    gr = purrr::map(gr, IRanges::resize, window * 2, fix = "center"))

vcf_files <- here::here("results/2022_08_04_ukbiobank_finemapping",
  "vcf/out") %>%
  list.files(full.names = TRUE, pattern = ".vcf.gz") %>%
  stringr::str_subset(".csi", negate = TRUE) %>%
  stringr::str_subset(".tbi", negate = TRUE)

get_genotype <- function(gr, vcf_files, stats) {

  chrom <- GenomeInfoDb::seqnames(gr) %>%
    as.character()
  vcf_file <- vcf_files

  svp <- VariantAnnotation::ScanVcfParam(
    geno = "GP", which = gr)

  vcf <- VariantAnnotation::readVcf(vcf_file, "hg19", param = svp)
  # vcf <- VariantAnnotation::readVcf(vcf_file, "GRCh38", param = svp)

  gmat <- VariantAnnotation::genotypeToSnpMatrix(vcf, uncertain = TRUE)
  gmat <- gmat[["genotypes"]] %>%
    as("numeric")

  ## there are some significant snps outside the 250 kbps range around the
  ## top snp
  idx <- is.na(gmat)
  idx <- names(which(colSums(idx) == 0))

  gmat[, idx]

}

library(VariantAnnotation)
library(Matrix)
library(snpStats)
library(GenomicRanges)
library(ggbio)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(ggrepel)


# get significant snps
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- org.Hs.eg.db 
ensdb <- EnsDb.Hsapiens.v75


locus_zoom_plot <- function(locus, index, db, stats, vcf_files, mart) {

  # browser()
  # get snps in locus

  stats_gr <- with(stats,
    GenomicRanges::GRanges(seqnames = str_c("0", chrom),
      ranges = IRanges::IRanges(start = pos,  width = 1)))

  ov_query <- IRanges::findOverlaps(locus, stats_gr)

  stats <- stats[S4Vectors::subjectHits(ov_query), ]

  geno_mat <- get_genotype(locus, vcf_files, stats)
  ld_mat <- cor(geno_mat)

  stats %<>%
    dplyr::filter(id %in% rownames(ld_mat))

  ld_vec <- ld_mat[stats$id, "rs4588"]

  stats %<>%
    dplyr::mutate(ld = ld_vec)

  snp_labels <- stats %>%
    dplyr::filter(ld^2 >= .3 & p <= 5e-8) %>%
    dplyr::mutate(
      coord = stringr::str_c(chrom, pos, pos, sep = ":"))

  # bm_query <- biomaRt::getBM(
  #   attributes = c("refsnp_id", "chr_name", "chrom_start", "allele"),
  #   filters = c("chromosomal_region"),
  #   values = unique(snp_labels$coord),
  #   mart = mart) %>%
  #   dplyr::mutate(
  #     coord = stringr::str_c(chr_name, chrom_start, chrom_start,
  #       sep = ":")) %>%
  #   dplyr::filter(stringr::str_detect(refsnp_id, regex("^rs"))) %>%
  #   dplyr::select(refsnp_id, coord)

  snp_labels %<>%
    dplyr::mutate(plot_id = id)

  plocus <- GenomicRanges::GRanges(
    seqnames = "4", ranges = ranges(locus))

  gene_annot <- ggbio::autoplot(db,
      AnnotationFilter::GRangesFilter(plocus),
      names.expr = "gene_name", label.size = 4) +
      theme_bw() + theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 14),
        legend.key.width = unit(.05, "npc"),
        legend.position = "top")

  main_plot <- stats %>%
    ggplot(aes(x = pos, y = -log10(p), fill = ld^2)) +
      geom_point(aes(shape = index), size = 5, colour = "black") +
      scale_shape_manual(values = c(21, 23), guide = FALSE) +
      labs(fill = expression(r^2), y = expression(-log[10](p.value))) +
      geom_hline(yintercept = 8, linetype = 2) +
      scale_fill_distiller(palette = "RdYlGn",
        labels = scales::percent_format(1)) +
      geom_text_repel(data = snp_labels, aes(label = plot_id)) +
      theme_bw()

  main_plot <- main_plot +
    theme(
      legend.key.width = unit(.15, "npc"),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14))

  tracks(main_plot, gene_annot, heights = c(2.5, 1), xlim = plocus,
    theme = theme_bw() + theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.key.width = unit(.05, "npc"),
        legend.position = "top")) +
    scale_x_sequnit("kb")

}



pngdir <- here::here("results", "2022_08_04_ukbiobank_finemapping",
  "figs")
fs::dir_create(pngdir)

mart <- useEnsembl("ENSEMBL_MART_SNP","hsapiens_snp",
  "https://feb2014.archive.ensembl.org")


top_snp %<>%
  dplyr::mutate(
    lplot = purrr::map2(gr, id2, locus_zoom_plot, ensdb, stats0,
      vcf_files, mart))

top_snp %<>%
  dplyr::mutate(
    filename = glue::glue("{dr}/ukb_afroc_25hvd_{id2}_locuszoom.pdf",
      dr = pngdir),
    walk2(filename, lplot, ggsave,
      width = 8, height = 6, units = "in"))


