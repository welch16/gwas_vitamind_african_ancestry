pacman::p_load(magrittr, tidyverse, gwasrapidd, snpStats, broom)

sccs_vitd_rsids <-
  here::here("results/2020_10_20_vitd_fine_mapping",
    "qs/sccs_vitdbp_biomart_queries.qs") %>%
  qs::qread()
sccs_vitd_colsummary <-
  here::here(
    "results/2020_10_20_vitd_fine_mapping/qs/sccs_vitdbp_col_summary.qs") %>%
  qs::qread()

# want tables with columns:
# rsid gene chrom pos ref/alt eaf beta std-err p.value

## 1. generate data_frame vitd, age, bmi, sex, pc_1,...,pc_20, table1 snps

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


datasets <- sccs_vitd_rsids %>%
  rename(id = term) %>%
  mutate(id = str_replace_all(id, "\\_", ":")) %>%
  dplyr::distinct(chrom, id) %>%
  nest(ids = c(id))

datasets %<>%
  mutate(
    plink_files = map(chrom, ~ get_plink_files(files, .)),
    snpmat = map2(plink_files, ids,
     ~ snpStats::read.plink(
        bed = str_subset(.x, ".bed"),
        bim = str_subset(.x, ".bim"),
        fam = str_subset(.x, ".fam"),
        select.snps = .y$id)))

# need to get response, pca, and covariates
workdir <- here::here("results/2020_09_01_merge_impute_sccs/plink2")
files <- list.files(workdir, full.names = TRUE, pattern = "vitdbp") %>%
  str_subset("rsq50")

build_covariate <- function(plink_files, sccs_files) {

  fam <- plink_files %>%
    str_subset(".fam") %>%
    readr::read_tsv(col_names = FALSE)

  pheno <- sccs_files %>%
    str_subset("_log_pheno") %>%
    read_delim(" ") %>%
    select(-ends_with("FID"))

  covs <- sccs_files %>%
    str_subset("bmi_age") %>%
    str_subset("log", negate = TRUE) %>%
    read_delim(" ") %>%
    select(-ends_with("FID"))

  fam %>%
    rename(IID = X2, sex = X5) %>%
    select(IID, sex) %>%
    inner_join(pheno, by = "IID") %>%
    inner_join(covs, by = "IID") %>%
    select(IID, vitd, everything()) %>%
    mutate(sex = if_else(sex == 1, "m", "f"))

}



build_data <- function(snpmat, data_f) {

  genomat <- snpmat$genotype %>%
    as("numeric")

  genomat %<>%
    as.data.frame() %>%
    as_tibble(rownames = "IID")

  modeldata <- data_f %>%
    inner_join(genomat, by = "IID") %>%
    select(-IID) %>%
    select(vitd, everything()) %>%
    mutate(sex = if_else(sex == "m", 1, 0)) %>%
    na.omit()

  return(modeldata)
}

## 2. perform linear models
datasets %<>%
  mutate(
    sccs_files = map(chrom, ~ get_plink_files(files, .)),
    data_f = map2(plink_files, sccs_files, build_covariate),
    data_f = map2(snpmat, data_f, build_data),
    model = map(data_f, ~ lm(vitd ~ ., data = .x)),
    model = map(model, broom::tidy))

vitdbp_data <- datasets$data_f[[1]] %>%
  dplyr::select(vitd) %>%
  dplyr::summarize(
    log_mean = mean(vitd),
    log_sd = sd(vitd),
    log_median = median(vitd),
    log_mad = mad(vitd),
    mean = mean(exp(vitd)),
    sd = sd(exp(vitd)),
    median = median(exp(vitd)),
    mad = mad(exp(vitd)))
    
fs::dir_create(here::here("results", "zz_methods", "csv", "summary"))
    
vitdbp_data %>%
  readr::write_csv(here::here("results", "zz_methods", "csv", "summary",
    "vitdbp_distribution.csv"))

datasets$data_f[[1]] %>%
  ggplot(aes(vitd)) + geom_histogram(bins = 31)
