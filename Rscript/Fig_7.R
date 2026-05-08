library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyr)
library(dplyr)

# load PGS data for long_df

#######################
## UKBB
#######################


p_UKBB <- ggplot(long_df, aes(x = Sex, y = OR, color = Model)) +
  geom_point(size = 5, position = position_dodge2(width = 0.4, reverse = TRUE)) +  # Plot OR as points
  geom_errorbar(aes(ymin = lw, ymax = up), width = 0.4, position = position_dodge2(width = 0.4, reverse = TRUE)) +  # Add error bars
  scale_color_manual(values = c("#fccccb", "#b0d992")) +
  labs(y = y_labs) +
  coord_flip() +
  ggh4x::facet_grid2(Phenotype ~ dataset, scales = "free_x", independent = "x") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size=10),
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size=17),
        text = element_text(size=14),
        legend.position = "bottom")

long_df_UKBB <- long_df

#########################
## external
#########################


p_external <- ggplot(long_df, aes(x = Sex, y = OR, fill = Model, color = Model)) +
  geom_boxplot(width = 0.6, position = position_dodge2(width = 0.4, reverse = TRUE)) +  # Adjust box width for better separation
  ggh4x::facet_grid2(Phenotype ~ dataset, scales = "free_x", independent = "x") +
  coord_flip() +  # Flip x and y axes
  scale_x_discrete(labels = c("All", "Female", "Male")) +  # Simplify x-axis labels
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # Increase space between groups
  scale_fill_manual(values = alpha(c("#fccccb", "#b0d992"), 0.6)) +
  scale_color_manual(values = c("#fccccb", "#b0d992")) +
  labs(y = y_labs) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size=10),
        legend.title = element_blank(),
        strip.text = element_text(size=17),
        text = element_text(size=14),
        legend.position = "bottom"
  )

