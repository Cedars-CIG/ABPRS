library(dplyr)
library(data.table)
library(ggplot2)
library(ggokabeito)
library(stringr)
library(patchwork)

point_and_line_plot <- function(data,
                                title = NULL,
                                show_range = TRUE,
                                legend_position = "none",
                                y_lab = "predicted mean value") {
    data |>
        ggplot(aes(x = allele, y = median, group = cluster, color = cluster)) +
        geom_line() +
        {
            if (show_range) {
                geom_segment(
                    aes(y = lo, yend = hi, x = allele, xend = allele),
                    linewidth = 0.75,
                    alpha = 0.5
                )
            }
        } +
        geom_segment(
            aes(y = iqr_lo, yend = iqr_hi, x = allele, xend = allele),
            linewidth = 0.75,
            alpha = 0.5
        ) +
        geom_point(size = 3) +
        labs(
            caption = title,
            y = y_lab
        ) +
        scale_color_okabe_ito() +
        facet_grid(
            positive ~ cluster,
            scales = "free_y",
            drop = FALSE,
            labeller = labeller(
                positive = c(
                    "Positive" = "Positive Slope",
                    "Negative" = "Negative Slope"
                )
            )
        ) +
        theme_classic() +
        theme(
            plot.caption = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_blank(),
            legend.position = legend_position,
            panel.grid.major.y = element_line(color = "gray", linewidth = 0.5),
            strip.text.x = element_text(size = 9)
        )
}


pheno_list <- c("AD", "breast_cancer", "hypertension", "t2d", "bmi", "cholesterol", "hdl", "ldl")
plot_title_list <- c("AD", "Breast Cancer", "Hypertension", "T2D", "BMI", "Cholesterol", "HDL", "LDL")

p_list <- vector("list", length(pheno_list))

for (i in 1:length(pheno_list)){

    print(pheno_list[i])

    if (pheno_list[i] %in% c("bmi", "cholesterol", "hdl", "ldl")) {
        y_lab <- "Predicted Mean Value"
    } else {
        y_lab <- "Predicted Disease Probability"
    }
    df <- read.csv(paste0(args[1],pheno_list[i],args[2]),header=T)

    dt <- setDT(df)
    p_list[[i]] <- dt |>
        point_and_line_plot(title = plot_title_list[i], show_range = show_range, y_lab = y_lab)

}

p <- p_list[[1]] + p_list[[2]] + p_list[[3]] + p_list[[4]] + 
  p_list[[5]] + p_list[[6]] + p_list[[7]] + p_list[[8]] + 
  plot_layout(ncol=2)

