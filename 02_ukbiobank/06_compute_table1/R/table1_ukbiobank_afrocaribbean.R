library(magrittr)
library(tidyverse)

source(here::here("src/gwas_funs_redo.R"))

# load summary statistics

work_dir <- here::here("results/2021_07_01_ukbiobank_stratified_gwas",
  "plink2/stats")
work_files <- list.files(work_dir, full.names = TRUE) %>%
  stringr::str_subset("chr4") %>%
  stringr::str_subset("afroc")

eaf <- work_files %>%
  stringr::str_subset("afreq") %>%
  vroom::vroom() %>%
  dplyr::rename_with(snakecase::to_snake_case)

stats <- work_files %>%
  stringr::str_subset("glm") %>%
  vroom::vroom() %>%
  dplyr::rename_with(snakecase::to_snake_case)

rsid_snp_keys <- tibble::tribble(
  ~ rsid, ~ id,
  "rs7041", "4:72618334:A:C",
  "rs842998", "4:72627797:G:C",
  "rs842873", "4:72696465:C:G",
  "rs11731496", "4:72819988:G:T",
  "rs705117", "4:72608115:C:T")

eaf %<>%
  dplyr::filter(id %in% rsid_snp_keys$rsid)

stats %<>%
  dplyr::filter(id %in% rsid_snp_keys$rsid)


eaf %<>%
  dplyr::select(id, alt_freqs) %>%
  dplyr::rename(eaf = alt_freqs)

stats %<>%
  dplyr::select(chrom, id, pos, beta, se, obs_ct, p)

snps <- dplyr::inner_join(
  stats, eaf, by = "id")

# Bioconductor part
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)
library(Matrix)

window <- 10


snps %<>%
  tidyr::nest(data = c(chrom, pos)) %>%
  dplyr::mutate(
    gr = purrr::map(data, create_gr),
    gr = purrr::map(gr, IRanges::resize, window * 2, fix = "center"))

vcf_files <- here::here("results/2022_08_04_ukbiobank_finemapping",
  "vcf/out") %>%
  list.files(full.names = TRUE, pattern = "_range.vcf.gz") %>%
  stringr::str_subset(".csi", negate = TRUE) %>%
  stringr::str_subset(".tbi", negate = TRUE)


get_genotype <- function(gr, vcf_files, stats) {

  chrom <- GenomeInfoDb::seqnames(gr) %>%
    as.character()
  vcf_file <- vcf_files

  GenomeInfoDb::seqlevels(gr) <- "04"
  svp <- VariantAnnotation::ScanVcfParam(
    geno = "GP", which = gr)

  vcf <- VariantAnnotation::readVcf(vcf_file, "hg19", param = svp)
  rnms <- rownames(vcf)
  rnms_rsids <- stringr::str_split(rnms, ",") %>%
    purrr::map_chr(1)

  stats %<>%
    dplyr::mutate(vcf_id = purrr::map(id, ~ which(. == rnms_rsids))) %>%
    dplyr::filter(purrr::map_dbl(vcf_id, length) > 0) %>%
    dplyr::mutate(rnm = purrr::map_chr(vcf_id, ~ rnms[.[1]]))

  vcf0 <- vcf[stats$rnm, ]
  geno_mat <- geno(vcf0)[["GP"]]
  geno_mat <- as.list(geno_mat)
  out_mat <- purrr::map(geno_mat, which.max)

  attr(out_mat, "dim") <- attr(geno_mat, "dim")
  attr(out_mat, "dimnames") <- attr(geno_mat, "dimnames")

  out_mat <- matrix(unlist(out_mat), ncol = ncol(vcf0), byrow = FALSE)
  rownames(out_mat) <- rownames(vcf0)
  colnames(out_mat) <- colnames(vcf0)


  ## there are some significant snps outside the 250 kbps range around the
  ## top snp
  out_mat <- out_mat - 1
  t(out_mat[stats$rnm, ])

}


snps %<>%
  dplyr::mutate(
    geno = purrr::map(gr, safely(get_genotype), vcf_files, snps))

geno <- snps %>%
  dplyr::filter(purrr::map_lgl(geno, ~ ! is.null(.$result))) %>%
  dplyr::mutate(
    geno = purrr::map(geno, "result"))

# add covariates and response
folder <- "2021_07_01_ukbiobank_stratified_gwas"
data_dir <- here::here("results", folder, "qs")
uk_data <- file.path(data_dir, "ukbiobank_grant_vars.qs") %>%
  qs::qread()

# uk_data %>%
#   dplyr::filter(pop == "afrocaribbean") %>%
#   dplyr::pull(vitamin_d) %>%
#   median(na.rm = TRUE)

# uk_data %>%
#   dplyr::filter(pop == "afrocaribbean") %>%
#   dplyr::pull(vitamin_d) %>%
#   {
#     quantile(., .75, na.rm = TRUE) -
#     quantile(., .25, na.rm = TRUE)}



# median(vitd_data$vdbp, na.rm = TRUE)
# iqr(vitd_data$vdbp, na.rm = TRUE)



princomp <- file.path(data_dir, "ukbiobank_princomp.qs") %>%
  qs::qread()
  
princomp %<>%
  dplyr::mutate(
    eid = stringr::str_split(iid, "\\_") %>%
      purrr::map_chr(1) %>%
      as.integer()) %>%
  dplyr::select(-iid)
  
uk_data %<>%
  dplyr::inner_join(princomp, by = "eid")

all_data <- uk_data %>%
  dplyr::filter(pop == "afrocaribbean") %>%
  dplyr::select(eid, vitamin_d, age_bin, sex, bmi, starts_with("pc"))

vec_to_tibble <- function(vec, name) {

  namevar <- rlang::sym(name)

  if (is.matrix(vec)) {

    nms <- colnames(vec)
    vec <- as.vector(vec)
    names(vec) <- nms

  }

  tibble::tibble(iid = names(vec), x = vec) %>%
    dplyr::rename(!! namevar := x)

}

geno_tibb <- purrr::map2(geno$geno, geno$id, vec_to_tibble) %>%
  purrr::reduce(dplyr::inner_join, by = "iid") %>%
  dplyr::rename(eid = iid) %>%
  dplyr::mutate(eid = as.integer(eid))

all_data %<>%
  dplyr::inner_join(geno_tibb, by = "eid")

snps %>%
  dplyr::select(-geno, -gr) %>%
  tidyr::unnest(cols = c(data)) %>%
  dplyr::mutate(
    pve = pve(beta, se, eaf, obs_ct),
    pve = scales::scientific(pve)) %>%
  readr::write_tsv(
    here::here("results", "2022_11_21_extra_tables",
      "tsv", "ukbiobank_afrocaribbean_individual_snps.tsv"))

all_data %>%
  names()

data1 <- all_data %>%
  dplyr::select(-rs705117)
data2 <- all_data

snp_stats <- snps %>%
  dplyr::select(id, data, obs_ct, eaf) %>%
  tidyr::unnest(cols = c(data))

model1 <- lm(log1p(vitamin_d) ~ ., data = data1) %>%
  broom::tidy() %>%
  dplyr::filter(str_detect(term, "rs")) %>%
  dplyr::rename(id = term) %>%
  dplyr::inner_join(snp_stats, by = "id") %>%
  dplyr::mutate(
    pve = pve(estimate, std.error, eaf, obs_ct),
    pve = scales::scientific(pve))

model2 <- lm(log1p(vitamin_d) ~ ., data = data2) %>%
  broom::tidy() %>%
  dplyr::filter(str_detect(term, "rs")) %>%
  dplyr::rename(id = term) %>%
  dplyr::inner_join(snp_stats, by = "id") %>%
  dplyr::mutate(
    pve = pve(estimate, std.error, eaf, obs_ct),
    pve = scales::scientific(pve))
