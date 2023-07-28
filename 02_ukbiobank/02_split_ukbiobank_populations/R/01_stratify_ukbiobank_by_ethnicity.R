pacman::p_load(
  magrittr,
  tidyverse,
  qs,
  rlang,
  Polychrome)

colo_cohort <- qs::qread(
  here::here("results",
    "2020_06_08_gwas_analysis",
    "qs",
    "ukbiobank_colorectal_cancer_2020_03_02.qs"))

cases_meaning <- qs::qread(
  here::here("results",
  "2020_06_08_gwas_analysis",
  "qs",
  "cases_meaning.qs"))

pca_results <- here::here("results",
  "2020_05_29_principal_component_analysis", "qs",
  "pca_results_plink2.qs") %>% qs::qread()

pca_results[["vectors"]] %<>%
  inner_join(colo_cohort, by = c(id = "eid")) %>%
  mutate(ethnicity = fct_explicit_na(ethnicity, "not-available"))

pca_results[["values"]] %<>%
  rename(eig_values = X1) %>%
  mutate(id = row_number())

pca_results[["vectors"]] %<>%
  dplyr::filter(ethnicity != "missing") %>%
  mutate(
    across(starts_with("pc"), list(~ . / sd(., na.rm = FALSE)),
      .names = "{.col}"))

comps <- pca_results[["vectors"]]

save_file <- function(data, filename) {

  data %>%
    dplyr::select(id) %>%
    dplyr::mutate(fid = id) %>%
    rlang::set_names("#FID", "IID") %>%
    readr::write_tsv(filename)

}


## european
folder <- "2021_01_06_ukbiobank_split_populations"
comps %>%
  dplyr::filter(
    ethnicity %in% c("british", "irish", "white",
    "any other white background")) %>%
    save_file(here::here("results", folder, "plink2",
      "european_samples.tsv"))

## african
comps %>%
  dplyr::filter(
    ethnicity %in% c("african")) %>%
    save_file(here::here("results", folder, "plink2",
      "african_samples.tsv"))

## caribbean
comps %>%
  dplyr::filter(
    ethnicity %in% c("caribbean")) %>%
    save_file(here::here("results", folder, "plink2",
      "caribbean_samples.tsv"))

## afro-caribbean
comps %>%
  dplyr::filter(
    ethnicity %in% c("african", "caribbean",
      "any other black background")) %>%
    save_file(here::here("results", folder, "plink2",
      "afrocaribbean_samples.tsv"))

## chinese and south-asian
asian_eth <- c("any other asian background", "bangladeshi", "chinese",
  "indian", "pakistani")

pca_mat <- comps %>%
  dplyr::filter(ethnicity %in% asian_eth) %>%
  dplyr::select(id, starts_with("pc")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("id") %>%
  as.matrix()

dist <- dist(pca_mat)
hc <- hclust(dist)

clusters <- cutree(hc, k = 3) %>%
  tibble(id = names(.), cluster = .)

comps %<>%
  dplyr::select(id, ethnicity) %>%
  dplyr::inner_join(clusters, by = "id") %>%
  dplyr::mutate(ethnicity = as.character(ethnicity))


comps %>%
  dplyr::select(-id) %>%
  table()

#                             cluster
# ethnicity                       1    2    3
#   any other asian background  306  863  276
#   bangladeshi                   0  190    2
#   chinese                    1238    1   14
#   indian                        0 4262   80
#   pakistani                     0 1310   14

comps %<>%
  dplyr::mutate(
    chinese = cluster == 1 | ethnicity == "chinese",
    south_asian = cluster %in% c(2, 3) &
      ethnicity != "chinese")

comps %>%
  dplyr::filter(chinese) %>%
    save_file(here::here("results", folder, "plink2",
      "chinese_samples.tsv"))

comps %>%
  dplyr::filter(south_asian) %>%
    save_file(here::here("results", folder, "plink2",
      "south_asian_samples.tsv"))

comps %>%
  save_file(here::here("results", folder, "plink2",
    "asian_samples.tsv"))
