library(magrittr)
library(tidyverse)
library(UWCCC.GWAS)

folder <- "2021_01_06_ukbiobank_split_populations"
data_dir <- here::here("results", folder, "qs")
out_dir <- here::here("results", folder, "figs")

qs_files <- list.files(data_dir, pattern = "qs")

make_m_plot <- function(data, outfile) {

  min_pval <- min(data$p)

  plot <- data %>%
    UWCCC.GWAS::manhattan_plot(sig = 5e-8, sig_color = "red") +
    labs(
      x = NULL,
      y = expression(-log[10](p))) +
    theme_minimal() +
    ylim(0, - 1.3 * log10(min_pval)) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))

  ggplot2::ggsave(
    filename = outfile,
    plot = plot,
    width = 12, height = 4.5, units = "in")

}

make_qplot <- function(data, outfile) {

  plot <- data %>%
    UWCCC.GWAS::qqplot_pvalue() +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_blank())

  ggplot2::ggsave(
    filename = outfile,
    plot = plot,
    width = 6, height = 6, units = "in")

}

qs_data <- purrr::map(file.path(data_dir, qs_files), qs::qread)

png_files <- stringr::str_replace(qs_files, ".qs", ".png")
png_files1 <- stringr::str_replace(png_files, ".png", "_manhattan.png")
png_files2 <- stringr::str_replace(png_files, ".png", "_qqplot.png")

map2(qs_data, file.path(out_dir, png_files1), make_m_plot)
map2(qs_data, file.path(out_dir, png_files2), make_qplot)
