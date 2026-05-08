library(ggplot2)
library(latex2exp)
library(readxl)

df <- read_xlsx("Fig_2b.xlsx")

labs <- c("Additive Genetic Model",
          "Dominant Genetic Model",
          "Recessive Genetic Model",
          "Co-dominant Genetic Model")
names(labs) <- c("add", "dom", "rec", "het")

df$coding <- factor(df$coding, levels=c("add", "dom", "rec", "het"))

ggplot(df, aes(x = x, y = Value, color=Model, fill=Model)) +
  geom_boxplot() +
  scale_x_discrete(labels = function(x) str_sub(x, 1, 2)) +
  facet_wrap(~coding, scales = "free",
             labeller = labeller(coding=labs)) + 
  labs(x="SNP Alleles", y="Predicted Disease Probability") + 
  scale_color_manual(labels = c("Additive Coding", TeX("${\\theta}_{SNP}$ Coding"), "Recaptured Coding"), 
                     values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  scale_fill_manual(labels = c("Additive Coding", TeX("${\\theta}_{SNP}$ Coding"), "Recaptured Coding"),
                    values = alpha(c("#66c2a5", "#fc8d62", "#8da0cb"), 0.6)) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
