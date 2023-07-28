
pacman::p_load(magrittr, tidyverse, gwasrapidd, snpStats, broom)

prev_studies <-
  here::here("results/zz_methods/qs/previous_studies_results.qs") %>%
  qs::qread()
sccs_vitd_annots <-
  here::here("results/zz_methods/qs/sccs_vdbp_significant_with_genes.qs") %>%
  qs::qread()
sccs_vitd_finemap <-
  here::here("results/2020_10_20_vitd_fine_mapping/qs",
    "fwd_step_reg_fine_mapping.qs") %>%
  qs::qread()
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

# piece 1: our data
base <- sccs_vitd_rsids %>%
  dplyr::select(-rsid, -alleles) %>%
  tidyr::unnest(cols = c(biomart_query)) %>%
  rename(id = term)

gene_annots <- sccs_vitd_annots %>%
  dplyr::select(id, genes) %>%
  tidyr::unnest(cols = c(genes)) %>%
  mutate(id = str_replace_all(id, ":", "\\_"))

sccs_vitd_colsummary %<>%
  dplyr::select(id, RAF, MAF) %>%
  mutate(id = str_replace_all(id, ":", "\\_"))

our_table <- base %>%
  inner_join(gene_annots, by = "id") %>%
  inner_join(sccs_vitd_colsummary, by = "id") %>%
  filter(map_int(str_split(allele, "\\/"), length) < 3) %>%
  filter(!str_detect(allele, "-")) %>%
  dplyr::select(refsnp_id, symbol, chrom, pos, allele, MAF, estimate,
    std.error, p.value) %>%
  dplyr::rename(
    rsid = refsnp_id,
    gene = symbol,
    beta = estimate,
    `ref/alt` = allele)

# piece 2: old data
prev_studies %<>%
  filter(author == "moy")

assocs <- prev_studies %>%
  pluck("assoc", 1)
variants <- prev_studies %>%
  pluck("variants", 1)

assocs %>% slotNames()

moy_study <- list(
  assocs %>% pluck("associations"),
  assocs %>% pluck("loci"),
  assocs %>% pluck("risk_alleles"),
  assocs %>% pluck("genes")) %>%
  reduce(partial(inner_join, by = "association_id")) %>%
  inner_join(variants@variants, by = "variant_id")

# rsid gene chrom pos ref/alt eaf beta std-err p.value
moy_study %<>%
  select(variant_id, gene_name, chromosome_name, chromosome_position,
    risk_allele, risk_frequency, beta_number, standard_error, pvalue) %>%
  rename(
    rsid = variant_id, gene = gene_name, chrom = chromosome_name,
    pos = chromosome_position, MAF = risk_frequency, beta = beta_number,
    std.error = standard_error, p.value = pvalue) %>%
  mutate(
    MAF = c(0.35, 0.44, 0.13),
    effect_allele = c("A", "C", "C"),
    risk_allele = str_c(effect_allele, risk_allele, sep = "/")) %>%
  rename(`ref/alt` = risk_allele) %>%
  add_row(
    rsid = "rs668443",
    gene = "GALNT2",
    chrom = "1",
    pos = 228276413,
    `ref/alt` = "A/T",
    MAF = 0.16,
    beta = 457.73,
    std.error = 109.81,
    p.value = 1.36e-5)


# piece 3: joint model

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

## 3. format table
datasets %<>%
  mutate(
    model = map(model, rename, id = term),
    model = map(model, filter, !str_detect(id, "pc")),
    model = map(model, filter, !id %in% c("sex", "age", "bmi")),
    model = map(model, filter, id != "(Intercept)"),
    model = map(model, mutate, id = str_remove_all(id, "`")))

joint <- datasets %>%
  select(chrom, model) %>%
  unnest(cols = c(model)) %>%
  mutate(id = str_replace_all(id, ":", "\\_"))


joint %<>%
  inner_join(gene_annots, by = "id") %>%
  inner_join(sccs_vitd_colsummary, by = "id") %>%
  inner_join(
    select(base, -estimate, -statistic, -p.value, -std.error),
      by = c("chrom", "id")) %>%
  filter(map_int(str_split(allele, "\\/"), length) < 3) %>%
  filter(!str_detect(allele, "-")) %>%
  dplyr::select(refsnp_id, symbol, chrom, pos, allele, MAF, estimate,
    std.error, p.value) %>%
  dplyr::rename(
    rsid = refsnp_id,
    gene = symbol,
    beta = estimate,
    `ref/alt` = allele)


# save stuff

our_table %>%
  write_csv(
    here::here("results", "zz_methods", "csv", "sccs_vitd_associations.csv"))

moy_study %>%
  write_csv(
    here::here("results", "zz_methods", "csv", "moy_vitd_associations.csv"))

joint %>%
  write_csv(
    here::here("results", "zz_methods", "csv", "sccs_vitd_assoc_joint.csv"))

list(
  "sccs_vitd" = our_table,
  "sccs_vitd_joint" = joint,
  "moy_vitd" = moy_study) %>%
  qs::qsave(
    here::here("results", "zz_methods", "qs", "vitd_associations.qs"))
