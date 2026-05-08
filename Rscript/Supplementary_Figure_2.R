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
library(readxl)
#setwd("~/your file path")

# 方法名称
methods <- c("Pre-trained PRS (LDpred)", "Lasso", "AB-PRS", "AB-PRS (single val)", "BHq")
methods_fac <- c( "AB-PRS", "AB-PRS (single val)", "BHq", "Lasso","Pre-trained PRS (LDpred)")
baseline<-"Pre-trained PRS (LDpred)"
# 自定义颜色
custom_colors <- c(
  "Pre-trained PRS (LDpred)" = "#F4A7B9",
  "Lasso" = "#d95f02",
  "AB-PRS" = "#A8D08D",
  "AB-PRS (single val)" = "#e7298a",
  "BHq" = "#7570b3"
)



all_df <- read_excel("Supplementary_Figure_2.xlsx")
all_df$Method <- factor(all_df$Method, levels = methods_fac)

############################
# Paired difference tests vs baseline + build arrow annotations
############################
# wide format: one row per replicate within each facet
wide_df <- all_df %>%
  select(Heritability, Scenario, Rep, Method, Loss) %>%
  pivot_wider(names_from = Method, values_from = Loss)


ann_list <- list()

Hs <- unique(wide_df$Heritability)

for (h in Hs) {
  
  scenarios <- unique(wide_df$Scenario[wide_df$Heritability == h])
  
  for (sc in scenarios) {
    
    sub_df <- wide_df %>% filter(Heritability == h, Scenario == sc)
    
    for (m in setdiff(methods, baseline)) {
      
      diff <- sub_df[[m]] - sub_df[[baseline]]
      
      tt <- tryCatch(t.test(diff, mu = 0), error = function(e) NULL)
      if (is.null(tt)) next
      
      pval <- tt$p.value
      mean_diff <- mean(diff, na.rm = TRUE)
      
      label <- NA_character_
      
      if (!is.na(pval) && pval < 0.05) {
        
        if (pval < 0.001) stars <- "***"
        else if (pval < 0.01) stars <- "**"
        else stars <- "*"
        
        # format p-value, e.g. 1.2e-05
        p_str <- formatC(pval, format = "e", digits = 0)
        
        if (mean_diff > 0) {
          label <- paste0(p_str, "\n+", stars)
        } else if (mean_diff < 0) {
          label <- paste0(p_str, "\n-", stars)
        }
      }else{
        label <- paste0(round(pval, digits = 3))
      }
      
      
      if (!is.na(label)) {
        
        y_method_max <- max(sub_df[[m]], na.rm = TRUE)
        
        ann_list[[length(ann_list) + 1]] <- data.frame(
          Heritability = h,
          Scenario = sc,
          Method = m,
          y = y_method_max * 1.03,
          arrow = label,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

ann_df <- bind_rows(ann_list)
ann_df$Method <- factor(ann_df$Method, levels = methods_fac)

############################
# Plot (your plot + arrow layer)
############################
p <- ggplot(all_df, aes(x = Method, y = Loss, fill = Method, color = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = alpha(custom_colors, 0.6)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = NULL, y = expression(R^2)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.placement = "outside",
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_blank()
  ) +
  facet_grid(Heritability ~ Scenario, scales = "free_y") +
  # give extra top space so arrows won't be clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  # add arrows
  geom_text(
    data = ann_df,
    aes(x = Method, y = y, label = arrow),
    inherit.aes = FALSE,
    size = 2.5
  )

print(p)
