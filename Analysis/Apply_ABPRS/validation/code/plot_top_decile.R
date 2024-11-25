# libraries and such -----------------------------------------------------------
library(data.table)
library(tidyverse)
library(patchwork)
library(ggokabeito)

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

bin_phenos  <- c("AD", "hypertension", "t2d", "breast_cancer")

pmbb <- map_dfr(
    bin_phenos,
    \(p) {
        fread(paste0("AB_PRS_result/pmbb/", p, "_V3_strat_or.csv")) |>
            mutate(pheno = p)
    }
) |> mutate(cohort = "PMBB")
emerge <- map_dfr(
    bin_phenos,
    \(p) {
        fread(paste0("AB_PRS_result/emerge/", p, "_V3_strat_or.csv")) |>
            mutate(pheno = p)
    }
) |> mutate(cohort = "eMERGE")

comb_plot_top0.1_bottom0.9 <- bind_rows(pmbb, emerge) |>
    mutate(model = factor(model, c("Pre-trained PRS", "AB-PRS")), cohort = factor(cohort, c("PMBB", "eMERGE"))) |>
    filter(range == "top0.1_mid0.0_0.9") |>
    ggplot(aes(x = sex, y = or, fill = model)) +
    labs(
        x = "", y = "Odds ratio", title = "Stratified PRS odds ratio by approach by outcome and sex",
        subtitle = "Top 10% vs. bottom 90%"
    ) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(cohort ~ pheno, scales = "free") +
    scale_fill_okabe_ito() +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank())
ggsave(
    plot = comb_plot_top0.1_bottom0.9, filename = "plots/comb_top0.1_bottom0.9_strat_or_plot.pdf", width = 9, height = 5, device = cairo_pdf
)

comb_plot_top0.1_bottom0.1 <- bind_rows(pmbb, emerge) |>
    mutate(model = factor(model, c("Pre-trained PRS", "AB-PRS")), cohort = factor(cohort, c("PMBB", "eMERGE"))) |>
    filter(range == "top0.1_bottom0.1") |>
    ggplot(aes(x = sex, y = or, fill = model)) +
    labs(
        x = "", y = "Odds ratio", title = "Stratified PRS odds ratio by approach by outcome and sex",
        subtitle = "Top 10% vs. bottom 10%"
    ) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(cohort ~ pheno, scales = "free") +
    scale_fill_okabe_ito() +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank())
ggsave(
    plot = comb_plot_top0.1_bottom0.1, filename = "plots/comb_top0.1_bottom0.1_strat_or_plot.pdf", width = 9, height = 5, device = cairo_pdf
)

comb_plot_top0.1_mid0.2 <- bind_rows(pmbb, emerge) |>
    mutate(model = factor(model, c("Pre-trained PRS", "AB-PRS")), cohort = factor(cohort, c("PMBB", "eMERGE"))) |>
    filter(range == "top0.1_mid0.4_0.6") |>
    ggplot(aes(x = sex, y = or, fill = model)) +
    labs(
        x = "", y = "Odds ratio", title = "Stratified PRS odds ratio by approach by outcome and sex",
        subtitle = "Top decile vs. middle 40th-60th percentiles"
    ) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(cohort ~ pheno, scales = "free") +
    scale_fill_okabe_ito() +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank())
ggsave(
    plot = comb_plot_top0.1_mid0.2, filename = "plots/comb_top0.1_mid0.2_strat_or_plot.pdf", width = 9, height = 5, device = cairo_pdf
)