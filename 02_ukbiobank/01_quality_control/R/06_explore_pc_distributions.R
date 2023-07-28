
pacman::p_load(
  magrittr,
  tidyverse,
  ukbtools)

principal_components <- qs::qread(
  here::here("results", "2020_05_29_principal_component_analysis", "qs",
    "pca_results_plink2.qs")) %>%
  pluck("vectors") %>%
  rename(eid = id)

# colo_cohort <- qs::qread(
#   here::here("results", "2020_06_08_gwas_analysis",
#     "qs", "ukbiobank_colorectal_cancer_2020_03_02.qs"))
#   rename(eid = id)

# colo_cohort %<>%
#   dplyr::mutate(eid = as.character(eid), fid = eid) %>%
#   dplyr::left_join(principal_components, by = "eid") %>%
#   dplyr::select(eid, fid, everything())

# NOTE: There was an issue with the logistic regression fit,
# probably due to the dimensions of the principal components


pc_long <- principal_components %>%
  tidyr::pivot_longer(-eid,
    names_to = "principal", values_to = "components") %>%
  dplyr::mutate(
    principal = fct_reorder(principal,
      as.numeric(str_remove(principal, "pc_"))))

pc_bp <- pc_long %>%
  ggplot(aes(principal, components)) + geom_boxplot()

ggsave(
  filename = here::here("results", "2020_06_08_gwas_analysis",
    "pc_dist.png"),
  pc_bp,
  width = 12,
  height = 4,
  units = "in")

pc_summary <- pc_long %>%
  group_by(principal) %>%
  summarize(
    mean = mean(components),
    sd = sd(components)
  )

pc_bp <- pc_long %>%
  inner_join(pc_summary, by = "principal") %>%
  ggplot(aes(principal, components / sd)) + geom_boxplot()

ggsave(
  filename = here::here("results", "2020_06_08_gwas_analysis",
    "pc_dist2.png"),
  pc_bp,
  width = 12,
  height = 4,
  units = "in")

