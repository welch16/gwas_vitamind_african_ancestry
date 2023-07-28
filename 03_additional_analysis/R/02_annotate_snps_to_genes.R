pacman::p_load(here, qs, magrittr, tidyverse, vroom, UWCCC.GWAS)

ukbiobank_vdbp <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "ukbiobank_vitd_stats.qs")) %>%
  dplyr::mutate(
    p = 2 * pnorm(t_stat, lower.tail = FALSE))

sccs_vdbp <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs")) %>%
  dplyr::mutate(
    p = 2 * pnorm(t_stat, lower.tail = FALSE))

sccs_hydroxy <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "sccs_hydroxy_stats.qs")) %>%
  dplyr::mutate(
    p = 2 * pnorm(t_stat, lower.tail = FALSE))


# filter to have significant snps
pval_thr <- 1e-4
ukbiobank_vdbp %<>% dplyr::filter(p <= pval_thr)
sccs_vdbp %<>% dplyr::filter(p <= pval_thr)
sccs_hydroxy %<>% dplyr::filter(p <= pval_thr)


# bioconductor stuff
pacman::p_load(
  EnsDb.Hsapiens.v75, AnnotationDbi,
  org.Hs.eg.db,
  TxDb.Hsapiens.UCSC.hg19.knownGene,
  BSgenome.Hsapiens.UCSC.hg19,
  biomaRt)


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- org.Hs.eg.db 

all_genes <- GenomicFeatures::genes(txdb,
  single.strand.genes.only = FALSE)
all_genes <- unlist(all_genes)

sccs_vdbp %<>% UWCCC.GWAS::annotate_genes(all_genes)
sccs_hydroxy %<>% UWCCC.GWAS::annotate_genes(all_genes)
ukbiobank_vdbp %<>% UWCCC.GWAS::annotate_genes(all_genes)

sccs_vdbp %>%
  qs::qsave(here::here("results", "zz_methods", "qs",
    "sccs_vdbp_significant_with_genes.qs"))

ukbiobank_vdbp %>%
  qs::qsave(here::here("results", "zz_methods", "qs",
    "ukbiobank_vdbp_significant_with_genes.qs"))

sccs_hydroxy %>%
  qs::qsave(here::here("results", "zz_methods", "qs",
    "sccs_hydroxy25_significant_with_genes.qs"))
