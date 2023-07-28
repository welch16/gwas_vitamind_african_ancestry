library(magrittr)
library(tidyverse)

source(here::here("src/gwas_funs_redo.R"))

# load summary statistics
stats1 <- ukb_afroc_gwas <-
  here::here("results/2022_08_04_ukbiobank_finemapping/qs",
   "afroc_ukb_stats_filter.qs") %>%
  qs::qread()

# fine-mapping parameters
pval <- 5e-8
window <- 250e3

stats1 %<>%
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
  list.files(full.names = TRUE, pattern = "_range.vcf.gz") %>%
  stringr::str_subset(".csi", negate = TRUE) %>%
  stringr::str_subset(".tbi", negate = TRUE)


get_genotype <- function(gr, vcf_files, stats) {

  chrom <- GenomeInfoDb::seqnames(gr) %>%
    as.character()
  vcf_file <- vcf_files

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

top_snp %<>%
  dplyr::mutate(
    geno = purrr::map(gr, get_genotype, vcf_files, stats1))

# add covariates and response
folder <- "2021_07_01_ukbiobank_stratified_gwas"
data_dir <- here::here("results", folder, "qs")
uk_data <- file.path(data_dir, "ukbiobank_grant_vars.qs") %>%
  qs::qread()

uk_data %>%
  dplyr::filter(pop == "afrocaribbean") %>%
  dplyr::pull(vitamin_d) %>%
  median(na.rm = TRUE)

uk_data %>%
  dplyr::filter(pop == "afrocaribbean") %>%
  dplyr::pull(vitamin_d) %>%
  {
    quantile(., .75, na.rm = TRUE) -
    quantile(., .25, na.rm = TRUE)}



median(vitd_data$vdbp, na.rm = TRUE)
iqr(vitd_data$vdbp, na.rm = TRUE)



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

fwd_step_regression <- function(geno, modeldata) {


  fwd_vars <- modeldata %>%
    dplyr::select(-eid, -vitamin_d) %>%
    names()

  geno_tib <- geno %>%
    as.data.frame() %>%
    tibble::rownames_to_column("eid") %>%
    tibble::as_tibble()

  modeldata %<>%
    dplyr::mutate(eid = as.character(eid))

  modeldata %<>%
    dplyr::inner_join(geno_tib, by = "eid") %>%
    dplyr::rename_all(list(~ str_replace_all(., ":", "\\_")))
  stop_criteria <- FALSE
  all_vars <- modeldata %>%
    dplyr::select(-vitamin_d, -eid) %>%
    names()

  linear_model <- function(var, fwd_vars, modeldata) {

    my_data <- modeldata %>%
      dplyr::select(tidyselect::one_of(c("vitamin_d", fwd_vars, var)))
    model <- lm(log1p(vitamin_d) ~ ., data = my_data)
    broom::tidy(model) %>%
      dplyr::filter(term == glue::glue("{var}"))
  }

  i <- 1
  modeldata %<>%
    dplyr::select(-eid)
  while (length(fwd_vars) < length(all_vars) & !stop_criteria) {

    message("iteration:", i)
    rem_vars <- all_vars[!all_vars %in% fwd_vars]
    message("linear models")
    pvalues <- furrr::future_map_dfr(rem_vars, linear_model, fwd_vars,
      modeldata)

    # select variable
    lowest <- pvalues %>%
      dplyr::top_n(1, -p.value) %>%
      dplyr::pull(term) %>%
      stringr::str_remove_all("`")
    pvalues %<>%
      dplyr::filter(term != glue::glue("{x}", x = lowest))

    stop_criteria <- all(pvalues$p.value > 5e-8)
    sum(pvalues$p.value <= 5e-8)
    fwd_vars <- c(fwd_vars, lowest)
    message(lowest)
    i <- i + 1
  }

  finaldata <- modeldata %>%
    dplyr::select(tidyselect::one_of(c("vitamin_d", fwd_vars)))
  finalmodel <- lm(log(vitamin_d) ~ ., data = finaldata)
  out <- broom::tidy(finalmodel) %>%
    dplyr::mutate(
      term = if_else(stringr::str_detect(term, regex("^pc")),
        term, stringr::str_replace_all(term, "\\_", ":")))

  return(out)

}

all_data <- uk_data %>%
  dplyr::select(eid, vitamin_d, age_bin, sex, bmi, starts_with("pc"))

top_snp %<>%
  dplyr::mutate(
    fwd = purrr::map(geno, fwd_step_regression, all_data))

model_snps <- top_snp %>%
  dplyr::select(id2, fwd) %>%
  tidyr::unnest(cols = c(fwd)) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::filter(! term %in% names(all_data)) %>%
  dplyr::filter(!stringr::str_detect(term, "age")) %>%
  dplyr::filter(!stringr::str_detect(term, "sex")) %>%
  dplyr::mutate(
    term = stringr::str_remove_all(term, "`"),
    coord = stringr::str_split(term, ":") %>%
      purrr::map_chr(~ stringr::str_c(.[1], .[2], .[2], sep = ":")))

model_snps %>%
  dplyr::inner_join(stats1, by = c(term = "id")) %>%
  dplyr::mutate(
    pve = pve(estimate, std.error, alt_freqs, 6934)) %>%
  qs::qsave(
    here::here("results", "2022_08_04_ukbiobank_finemapping", "qs",
      "ukb_25hvitd_fwd_step_range.qs"))
