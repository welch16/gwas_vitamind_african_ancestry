
pacman::p_load(magrittr, tidyverse, gwasrapidd)

# adding manually two papers
# Hong et al

tables <- tibble::tribble(
  ~ author, ~ study_id,
  "hong", "GCST005782",
  "moy", "GCST002414",
  "jiang", "GCST005366")

tables %<>%
  mutate(
    study = map(study_id, gwasrapidd::get_studies),
    assoc =  map(study_id, gwasrapidd::get_associations),
    variants = map(study_id, gwasrapidd::get_variants))

tables %>%
  qs::qsave(here::here("results", "zz_methods", "qs",
    "previous_studies_results.qs"))
