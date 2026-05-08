library(dplyr)
library(readr)
library(data.table)
library(reshape2)
library(patchwork)
library(ggplot2)

p <- ggplot(df, aes(x = group, fill = coding)) + 
  geom_bar() +
  scale_fill_manual(values = c("Additive" = "#7FB4A6", "Dominant" = "#E3D5AE", "Recessive" = "#82531E", "Co-dominant" = "#B5823E")) +
  facet_wrap(~geno, scales = "free", ncol = 4, 
             labeller = labeller(geno=labs)) + 
  labs(y = "Number of Selected SNP") + 
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
