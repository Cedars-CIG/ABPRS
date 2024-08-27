# libraries --------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(ggtext)

n              <- 8000 # number of SNPs to generate
output_path    <- "~/Downloads/"
plot_file_name <- "manhattan.png"

# specifications ---------------------------------------------------------------
cols <- c(
    "ADD"    = "#D55E00",
    "DOM"    = "#0072B2",
    "REC"    = "#CC79A7",
    "OVRDOM" = "#009E73",
    "OVRREC" = "#E69F00"
)

shps <- c(
    "ADD"    = 3,
    "DOM"    = 15,
    "REC"    = 16,
    "OVRDOM" = 17,
    "OVRREC" = 18
)

# data -------------------------------------------------------------------------
generate_skewed_pvalues <- function(n, shape1 = 10, shape2 = 6, power = 9) {
  # n: Number of p-values to generate
  # shape1, shape2: Shape parameters for the Beta distribution (default gives a right-skewed distribution)
  # power: Exponent to which p-values are raised to create more extreme values (default is 1, no change)
  
  # Generate p-values from a Beta distribution
  p_values <- rbeta(n, shape1, shape2)
  
  # Apply the power transformation to make values more extreme
  p_values <- p_values^power
  
  # Ensure p-values are between 0 and 1
  p_values <- pmin(pmax(p_values, 0), 1)
  
  return(p_values)
}

generate_data <- function(n) {
    data.table(
        snp     = rep(paste0("rs", sample(100000:999999, size = n, replace = FALSE)), 5),
        model   = rep(c("ADD", "DOM", "REC", "OVRDOM", "OVRREC"), each = n),
        p_value = unlist(map(1:5, \(x) generate_skewed_pvalues(n)))
    )
}

process_data <- function(data, threshold = 0.05/n) {
    # Additive p-values
    add_data <- data[model == "ADD", ]
    alt_data <- data[model != "ADD", ]
    # Alternate (non-additive) p-values
    alt_data <- alt_data[alt_data[, .I[which.min(p_value)], snp][["V1"]], ]
    # stack
    stack <- rbindlist(list(add_data, alt_data))
    min_idx <- stack[, .I[which.min(p_value)], snp][["V1"]]
    stack[, `:=` (min = 0, id = as.numeric(factor(snp)), log10p = -log10(p_value))]
    stack[min_idx, min := 1]
    stack[, color := NA_character_]
    stack[min == 1 & p_value < threshold, color := model][]
}

color_man_plot <- function(
    data,
    threshold  = NULL,
    x_lab      = "SNP",
    y_lab      = "-log10(p-value)",
    point_size = 2,
    alpha      = 0.75
) {
    data |>
        ggplot(aes(x = id, y = log10p, color = color, shape = model)) +
        geom_point(size = point_size, alpha = alpha) +
        {if (!is.null(threshold)) geom_hline(yintercept = threshold, color = "red")} +
        labs(x = x_lab, y = y_lab) +
        xlim(0, data[, length(unique(id))]) +
        scale_color_manual(breaks = unique(data[, model]), values = cols, name = "model") +
        scale_shape_manual(name = "model", breaks = unique(data[, model]), values = shps) +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        )
}

# plot -------------------------------------------------------------------------
full_plot <- generate_data(8000) |>
    process_data(threshold = 0.05/8000) |>
    color_man_plot(threshold = -log10(0.05/8000)) +
    labs(
        caption = paste0("**NOTES**:<br>",
        " - n = 8,000 SNPs<br>",
        " - Two p-value per SNP, one for the additive model and one for the best non-additive model.<br>",
        " - Points are colored by the best model for each SNP, provided it meets the significance threshold.<br>",
        " - Red line indicates the Bonferroni-corrected significance threshold."
        )
    ) +
    theme(
        plot.caption = element_markdown(hjust = 0)
    )

# save plots -------------------------------------------------------------------
ggsave(
    plot     = full_plot,
    filename = paste0(output_path, plot_file_name),
    width    = 10,
    height   = 5.5,
    dpi      = 240,
    units    = "in"
)
