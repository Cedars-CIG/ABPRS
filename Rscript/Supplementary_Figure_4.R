library(tidyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggh4x)
library(writexl)
#setwd("~/your file path")

############################
##### pre-trained PRS LD ### 
############################

methods <- c("Pre-trained PRS (basic)", "Lasso", "AB-PRS", "AB-PRS (single val)", "BHq", "Bayesian")


custom_colors <- c(
  "Pre-trained PRS (basic)" = "#F4A7B9",
  "Lasso" = "#d95f02",
  "AB-PRS" = "#A8D08D",
  "AB-PRS (single val)" = "#e7298a",
  "BHq" = "#7570b3",
  "Bayesian"="#85A9CE"
)


all_df <- read_excel("Supplementary_Figure_4_left.xlsx")
#####################
# plot
#####################


p1 <- ggplot(all_df, aes(x = Method, y = Loss, fill = Method, color = Method)) + 
  geom_boxplot(show.legend = FALSE) + 
  geom_point(aes(fill = Method), size = -1) +
  scale_fill_manual(values = alpha(custom_colors, 0.4)) + 
  scale_color_manual(values = custom_colors) + 
  theme_minimal() + 
  labs(x = NULL, y = NULL) +
  guides(
    fill = guide_legend(override.aes = list(size = 3, shape = 15))
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.placement = "outside",
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_blank(),
    strip.text.y = element_blank()   # 直接去掉右边的 strip 文本
  ) +
  facet_grid(
    Heritability ~ Scenario,
    scales = "free_y",
    labeller = labeller(
      Scenario = function(x) "R²"   # 顶部 strip 全部换成 R²
    )
  ) 

legend_box <- ggpubr::get_legend(p1)

p1 <- p1 + theme(legend.position = "none") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x))

p1


##### selection propotion#######


custom_colors <- c(
  "Lasso" = "#d95f02",
  "AB-PRS" = "#A8D08D",
  "AB-PRS (single val)" = "#e7298a",
  "BHq" = "#7570b3",
  "Bayesian"="#85A9CE"
)


all_df <- read_excel("Supplementary_Figure_4_right.xlsx")
#####################
# plot
#####################
p2 <- ggplot(all_df, aes(x = group, y = value, fill=Method, color = Method)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = alpha(custom_colors, 0.4)) +
  scale_color_manual(values = custom_colors)+
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "grey70"),  # 行、列共用
    strip.placement = "outside",  # 让标题放在外侧（可选）
    strip.text = element_text(size = 11, face = "bold"),
    #panel.border = element_blank(), # 不要四周边框
    panel.grid = element_blank(),
    #axis.line.x = element_line(color = "black"), # 只保留下边框
    #axis.line.y = element_line(color = "black")  # 只保留左边框
    panel.border = element_rect(color = "black", fill = NA), # 每个子图加边框
    axis.line = element_blank() 
  ) +
  facet_grid(Heritability ~ Scenario, scales = "free_y")

p2<-p2+scale_x_discrete(labels = function(x) gsub(" ", "\n", x))
p2

layout <- "
AB
cc
"

p<- p1+p2 + legend_box+
  plot_layout(design = layout, heights = c(20, 1)
              ,widths = c(2, 2))
p
