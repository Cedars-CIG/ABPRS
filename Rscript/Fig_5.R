library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(latex2exp)
library(scales)
library(ggh4x)


# Define shapes and colors
shapes <- c("Additive" = 15, "Dominant" = 16, "Recessive" = 17, "Co-dominant" = 18)
colors <- c("Additive" = "#7FB4A6", "Dominant" = "#E3D5AE", "Recessive" = "#82531E", "Co-dominant" = "#B5823E", "Sub-optimal" = "grey")


# Custom transformation function
custom_trans <- function(y0) {
  trans_new(
    name = "custom",
    transform = function(y) ifelse(y <= y0, y, y0 + (y - y0) / 50),
    inverse = function(y) ifelse(y <= y0, y, y0 + (y - y0) * 50)
  )
}


# break point
y0 <- 20

# Plot the data
p <- ggplot(df_adjust, aes(x = BPcum, y = pval, shape = type, color = optimal_type, alpha = is_max)) +
  geom_point(size = 3) +
  facet_wrap(~pheno, nrow = 8, scales = "free", strip.position="right") + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(labels = c("Additive", "Dominant", "Recessive", "Co-dominant",  "Sub-optimal Additive"),
                     values = colors) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.8)) +  # Set alpha values based on is_max
  labs(x = "Chromosome", y = TeX("$-log_{10}(p-value)$"), color = "") +
  guides(shape = "none",  # Remove shape legend
         color = guide_legend(override.aes = list(shape = c(15, 16, 17, 18, 15))),
         alpha = "none") +  # Remove alpha from legend
  scale_x_continuous(label = df_bp_info$CHR, breaks= df_bp_info$BP_mid) +
  scale_y_continuous(trans = custom_trans(y0),
                     breaks = function(y) {
                      if (max(y) > 100) return(c(seq(0, y0-5, by = 5), seq(y0+80, max(y), by = 200))) else {return(c(seq(0, y0-5, by = 5), 70))}
                     }) + 
  geom_vline(xintercept=c(df_bp_info$BP_start, df_bp_info$BP_end[22]), linetype="dashed", colour="grey") + 
  annotate("text", y = c(y0-0.2, y0+20), x = -Inf, label = "~", size=unit(6,"pt")) + coord_cartesian(clip = "off") + 
  #geom_rect(aes(xmin=-Inf, xmax=+Inf, ymin=y0-0.4, ymax=y0+2, fill="white", color="white")) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        strip.text = element_text(angle=0, hjust=0.5),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text()
        ) + 
  guides(y = guide_axis_truncated(
          trunc_lower = c(-Inf, y0+20),
          trunc_upper = c(y0-0.2, Inf)
        ))
