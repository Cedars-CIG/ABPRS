# libraries --------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(ggokabeito)

data_path      <- "https://raw.githubusercontent.com/wuoli/ABPRS/main/Analysis_clustering/breast_cancer_PGS000004_kmeans_rslt.csv?token=GHSAT0AAAAAACWGUKS7BVPU7BZK5X5DJYR2ZWMZAEA" # update
plot_title     <- "Breast cancer" # update
show_range     <- TRUE
output_path    <- "~/Downloads/" # update
plot_file_name <- "breast_cancer_PGS000004_kmeans_rslt.png" # update

# load data --------------------------------------------------------------------
d <- fread(data_path)

# functions --------------------------------------------------------------------
## process data
process_data <- function(data) {
    data |>
        data.table::melt(
            id.vars = c("SNP", "cluster"),
            measure.vars = c("p_AA", "p_Aa", "p_aa"),
            variable.name = "allele",
            value.name = "value"
        ) |>
        (\(x) x[, .(
            mean = mean(value),
            median = median(value),
            iqr_lo = quantile(value, 0.25),
            iqr_hi = quantile(value, 0.75),
            lo = min(value),
            hi = max(value)
            ), c("allele", "cluster")])() |>
        dplyr::mutate(
            allele = factor(data.table::fcase(
                allele == "p_AA", 1,
                allele == "p_Aa", 2,
                allele == "p_aa", 3
            ), levels = c(1, 2, 3), labels = c("AA", "Aa", "aa")),
            cluster = as.factor(cluster)
        )
}

## plot data
point_and_line_plot <- function(data, title = NULL, show_range = TRUE, legend_position = "none") {
    data |>
        ggplot(aes(x = allele, y = median, group = cluster, color = cluster)) +
        geom_line() +
        {if (show_range) geom_segment(aes(y = lo, yend = hi), linewidth = 0.75, alpha = 0.5)} +
        geom_segment(aes(y = iqr_lo, yend = iqr_hi), linewidth = 1.75, alpha = 0.75) + # add IQR
        geom_point(size = 3) +
        labs(
            title = title,
            x = "",
            y = "Value"
        ) +
        scale_color_okabe_ito() +
        facet_wrap(~cluster, nrow = 1) + # modify to add facet for positive/negative
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0, face = "bold", size = 15),
            legend.position = legend_position,
            panel.grid.major.y = element_line(color = "gray", size = 0.5)
        )
}


# process data -----------------------------------------------------------------
example_plot <- d |>
    process_data() |>
    point_and_line_plot(title = plot_title, show_range = show_range)

ggsave(
    paste0(output_path, plot_file_name),
    plot = example_plot,
    width = 7,
    height = 4,
    dpi = 240,
    units = "in"
)
