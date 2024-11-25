# libraries and such -----------------------------------------------------------
library(data.table)
library(tidyverse)
library(patchwork)
library(ggokabeito)
library(DescTools)
library(patchwork)

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

bin_phenos  <- c("AD", "hypertension", "t2d", "breast_cancer")

pmbb <- map_dfr(
    bin_phenos,
    \(p) {
        fread(paste0("AB_PRS_pred/pmbb/", p, "_pred.txt")) |>
            mutate(pheno = p)
    }
) |> mutate(cohort = "PMBB")

emerge <- map_dfr(
    bin_phenos,
    \(p) {
        fread(paste0("AB_PRS_pred/emerge/", p, "_pred.txt")) |>
            mutate(pheno = p)
    }
) |> mutate(cohort = "eMERGE")

# functions --------------------------------------------------------------------
make_deciles <- function(x, p) {
    x |>
        filter(pheno == p) |>
        mutate(
            decile1 = ntile(pred1, 10),
            decile2 = ntile(pred2, 10),
        )
}
print_catt_p <- function(x) {
    tmp_p <- x$p.value
    if (tmp_p < 0.001) {
        "P < 0.001"
    } else {
        paste0("P = ", trimws(format(round(tmp_p, 3), nsmall = 3)))
    }
}


quick_prep <- function(x) {
    dt_p1 <- x |>
        summarize(
            case = sum(y == 1),
            control = sum(y == 0),
            .by = decile1
        ) |>
        mutate(
            pp = case / (case + control),
            p_se = qnorm(0.975) * sqrt((pp * (1-pp)) / (case + control)),
            p_lo = pp - p_se,
            p_hi = pp + p_se,
            decile = as.numeric(decile1),
            model = "Pre-trained PRS"
        ) |>
        select(
            decile, p = pp, p_lo, p_hi, model
        )
    dt_p2 <- x |>
        summarize(
            case = sum(y == 1),
            control = sum(y == 0),
            .by = decile2
        ) |>
        mutate(
            pp = case / (case + control),
            p_se = qnorm(0.975) * sqrt((pp * (1-pp)) / (case + control)),
            p_lo = pp - p_se,
            p_hi = pp + p_se,
            decile = as.numeric(decile2),
            model = "AB-PRS"
        ) |>
        select(
            decile, p = pp, p_lo, p_hi, model
        )
    bind_rows(
        dt_p1, dt_p2
    ) |> mutate(model = factor(model, c("Pre-trained PRS", "AB-PRS")))
}

quick_prep(pmbb_dec) |>
    pull(p_hi) |>
    max() |>
    (\(x) x* 1.05)()

quick_plot2 <- function(dt_data, .pheno, cohort) {

    dt <- dt_data |> filter(pheno == .pheno)
    # make dec
    dt_dec   <- make_deciles(dt, .pheno)

    # prep data
    plot_dat   <- quick_prep(dt_dec)
    n_cases    <- format(dt[, sum(y == 1)], big.mark = ",")
    n_controls <- format(dt[, sum(y == 0)], big.mark = ",")

    # calculate p
    dt_cat1 <- CochranArmitageTest(table(dt_dec[, decile1], dt_dec[, y]), alternative = "o")
    dt_cat2 <- CochranArmitageTest(table(dt_dec[, decile2], dt_dec[, y]), alternative = "o")
    dt_cat_p1 <- print_catt_p(dt_cat1)
    dt_cat_p2 <- print_catt_p(dt_cat2)
    dt_cat_p <- data.table(
        model = c("Pre-trained PRS", "AB-PRS"),
        p     = c(dt_cat_p1, dt_cat_p2),
        y     = c(
            plot_dat |>
                    pull(p_hi) |>
                    max() |>
                    (\(x) x* 1.15)(),
            plot_dat |>
                    pull(p_hi) |>
                    max() |>
                    (\(x) x* 1.05)()
        )
    ) |> mutate(model = factor(model, c("Pre-trained PRS", "AB-PRS")), cohort = cohort)

    # stack data
    plot_dat |>
        ggplot(aes(x = decile, y = p, fill = model)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(
            aes(ymin = p_lo, ymax = p_hi),
            width = 0.5, linewidth = 1, position = position_dodge(0.9),
        ) +
        geom_text(
            aes(x = 1, y = y, color = model, label = p),
            data = dt_cat_p, show.legend = FALSE
        ) +
        labs(
            title    = paste0("% ", .pheno, " cases by PRS decile in ", cohort),
            subtitle = paste0("n = ", n_cases, " cases; n = ", n_controls, " controls"),
            x        = "PRS decile",
            y        = "% cases",
            caption  = paste0("Notes:\n",
            "  - One-sided p-value for Cochran-Armitage test for trend reported")
        ) +
        scale_fill_okabe_ito() +
        scale_color_okabe_ito() +
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(breaks = 1:10) +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            plot.caption = element_text(hjust = 0),
            plot.title = element_text(face = "bold")
        )
}

patch_plots <- function(pmbb_dt, emerge_dt, pheno) {
    pmbb_plot   <- quick_plot2(dt_data = pmbb_dt, .pheno = pheno, cohort = "PMBB") +
        labs(title = "A. PMBB", caption = "") +
        theme(legend.position = "none")
    emerge_plot <- quick_plot2(dt_data = emerge_dt, .pheno = pheno, cohort = "PMBB") +
        labs(title = "B. eMERGE", caption = "") +
        theme(legend.position = "bottom")

    (pmbb_plot / emerge_plot) +
        plot_annotation(
            title = paste0("% ", pheno, " cases by PRS decile by cohort"),
            caption = paste0("Notes:\n",
            "  - One-sided p-value for Cochran-Armitage test for trend reported")) &
            theme(plot.caption = element_text(hjust = 0), plot.title = element_text(face = "bold"))
}

for (i in c("AD", "t2d", "breast_cancer", "hypertension")) {
    print(i)
    ggsave(
        plot = patch_plots(pmbb_dt = pmbb, emerge_dt = emerge, pheno = i),
        width = 7, height = 8, device = cairo_pdf,
        filename = paste0("plots/", i, "_decile_plot.pdf")
    )
}
