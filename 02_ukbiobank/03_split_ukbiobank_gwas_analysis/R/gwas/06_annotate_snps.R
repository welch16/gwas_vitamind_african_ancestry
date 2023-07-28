
library(here)
library(qs)
library(magrittr)
library(tidyverse)
# library(UWCCC.GWAS)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)


pval_thr <- 1e-4

folder <- "2021_07_01_ukbiobank_stratified_gwas"


asian_gwas <- qs::qread(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_asian.qs")) %>%
    filter(p <= pval_thr)

white_gwas <- qs::qread(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_white.qs")) %>%
    filter(p <= pval_thr)

afroc_gwas <- qs::qread(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_afrocaribbean.qs")) %>%
    filter(p <= pval_thr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- GenomicFeatures::genes(txdb,
  single.strand.genes.only = FALSE)
all_genes <- unlist(all_genes)

convert_gwas <- function(gwas) {

  gwas %>%
    dplyr::select(chr, pos, id) %>%
    dplyr::mutate(chrom = str_remove(chr, "chr"))

}

asian_gr <- convert_gwas(asian_gwas)
white_gr <- convert_gwas(white_gwas)
afroc_gr <- convert_gwas(afroc_gwas)

annotate_genes <- function(snps, all_genes, db) {

  gr <- with(snps,
    GenomicRanges::GRanges(
      seqnames = stringr::str_c("chr", chrom),
      ranges = IRanges::IRanges(start = pos, width = 1))) %>%
    magrittr::set_names(dplyr::pull(snps, id))

  nms <- names(gr)
  nearest_genes <- IRanges::nearest(
    gr, all_genes, select = "all", ignore.strand = TRUE)

  snp_ids <- S4Vectors::queryHits(nearest_genes)
  gene_ids <- S4Vectors::subjectHits(nearest_genes)
  entrezids <- names(all_genes[gene_ids])

  gene_search <- AnnotationDbi::select(db,
    columns = c("SYMBOL"),
    keytype = "ENTREZID", keys = unique(entrezids)) %>%
    tibble::as_tibble() %>%
    dplyr::rename_all(list(snakecase::to_snake_case))

  out <- tibble(id = nms[snp_ids], entrezid = entrezids) %>%
    dplyr::left_join(gene_search, by = "entrezid") %>%
    dplyr::select(-entrezid)

  id <- NULL
  db <- NULL
  entrezid <- NULL
  symbol <- NULL

  return(
    snps %>%
      dplyr::inner_join(out, by = "id") %>%
      tidyr::nest(genes = c(symbol)))
}

human <- org.Hs.eg.db

asian_gr %<>%
  annotate_genes(all_genes, human) %>%
  dplyr::select(-chrom, -chr, -pos)
white_gr %<>%
  annotate_genes(all_genes, human) %>%
  dplyr::select(-chrom, -chr, -pos)
afroc_gr %<>%
  annotate_genes(all_genes, human) %>%
  dplyr::select(-chrom, -chr, -pos)

asian_gr %>%
  qs::qsave(here::here("results", folder, "qs/ukb_asian_annot_snps.qs"))
white_gr %>%
  qs::qsave(here::here("results", folder, "qs/ukb_white_annot_snps.qs"))
afroc_gr %>%
  qs::qsave(here::here("results", folder, "qs/ukb_afroc_annot_snps.qs"))
