
library(magrittr)
library(tidyverse)

table_file <- here::here("results/2022_11_11_paper_review/data",
  "Table 1_UKB_strat_ancestry.xlsx")

snp_info <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "A4:F147") %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::select(snp)


revez <- list()
revez[["before_cojo"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "G4:I147") %>%
  dplyr::rename_with(snakecase::to_snake_case)
revez[["after_cojo"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "J4:L147") %>%
  dplyr::rename_with(snakecase::to_snake_case)
revez[["bmi_cond"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "M4:O147") %>%
  dplyr::rename_with(snakecase::to_snake_case)

stratified <- list()
stratified[["afrocaribbean"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "R4:AB147") %>%
  dplyr::rename_with(snakecase::to_snake_case)
stratified[["asian"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "AC4:AM147") %>%
  dplyr::rename_with(snakecase::to_snake_case)
stratified[["white"]] <- readxl::read_excel(
  table_file, sheet = "Table 1", range = "AN4:AX147") %>%
  dplyr::rename_with(snakecase::to_snake_case)

revez %<>%
  purrr::map(~ dplyr::bind_cols(snp_info, .))

stratified %<>%
  purrr::map(~ dplyr::bind_cols(snp_info, .))

make_scatter <- function(revez, ukb_strat) {

  revez_beta <- names(revez)[2]
  plot_data <- revez %>%
    dplyr::inner_join(ukb_strat, by = "snp") %>%
    dplyr::mutate(beta = as.numeric(beta)) %>%
    na.omit()

  cor_test <- cor.test(
    x = plot_data[[revez_beta]],
    y = plot_data[["beta"]], method = "spearman") %>%
    broom::tidy()

  annot <- glue::glue(
    "Spearman correlation: {rho}\n p.value = {pval}",
    rho = round(cor_test$estimate, 3),
    pval = scales::scientific(cor_test$p.value))

  annot_data <- data.frame(
    xpos = -Inf, ypos = Inf, hjustvar = 0, vjustvar = 1,
    annot = as.character(annot))

  plot_data %>%
    ggplot(aes_string(x = revez_beta, y = "beta")) +
      geom_point() +
      geom_smooth(se = FALSE, method = "lm", formula = y ~ x) +
      geom_text(
        data = annot_data,
        aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
          label = annot), size = 5) +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))
}


figsdr <- here::here("results/2022_11_11_paper_review/figs")

beta_values <- as.character(
  glue::glue("Revez: {rev} beta coefficients", rev = names(revez)))

plots <- purrr::map2(
  revez, beta_values,
  ~ make_scatter(.x, stratified[["afrocaribbean"]]) +
    labs(x = .y,
      y = "UK Biobank stratified Afrocaribbean beta coefficients"))

filenames <- file.path(
  figsdr,
  as.character(
    glue::glue(
      "scatterplot_{revez}_vs_ukb_strat_afrocaribbean.png",
      revez = names(revez))))

purrr::map2(
  filenames, plots,
  ggsave, width = 7, height = 7, units = "in", bg = "white")

plots <- purrr::map2(
  revez, beta_values,
  ~ make_scatter(.x, stratified[["asian"]]) +
    labs(x = .y,
      y = "UK Biobank stratified Asian beta coefficients"))

filenames <- file.path(
  figsdr,
  as.character(
    glue::glue(
      "scatterplot_{revez}_vs_ukb_strat_asian.png",
      revez = names(revez))))

purrr::map2(
  filenames, plots,
  ggsave, width = 7, height = 7, units = "in", bg = "white")

plots <- purrr::map2(
  revez, beta_values,
  ~ make_scatter(.x, stratified[["white"]]) +
    labs(x = .y,
      y = "UK Biobank stratified European beta coefficients"))

filenames <- file.path(
  figsdr,
  as.character(
    glue::glue(
      "scatterplot_{revez}_vs_ukb_strat_european.png",
      revez = names(revez))))

purrr::map2(
  filenames, plots,
  ggsave, width = 7, height = 7, units = "in", bg = "white")
