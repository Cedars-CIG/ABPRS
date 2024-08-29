library(ggplot2)
library(dplyr)
library(tidyr)

# Load Data
ad_cluster <- read.csv("/Users/oreomilk/ABPRS/Analysis_clustering/AD_DNAnexus_PGS001775_kmeans_rslt.csv")

# Step 1: Reorder the Clusters
# Wanling's Code

# Step 2: Distinguish Between Positive and Negative Slopes
df_summary_sep <- ad_cluster %>%
  mutate(slope_group = ifelse(slope1 > 0, "Positive Slope", "Negative Slope"))

# Step 3: Get Median, Quantiles, and Make Table Longer
df_summary_sep_iqr <- df_summary_sep %>%
  group_by(cluster, slope_group) %>%
  summarize(
    AA = median(p_AA, na.rm = TRUE),
    Q1_pAA = quantile(p_AA, 0.25, na.rm = TRUE),
    Q3_pAA = quantile(p_AA, 0.75, na.rm = TRUE),
    Aa = median(p_Aa, na.rm = TRUE),
    Q1_pAa = quantile(p_Aa, 0.25, na.rm = TRUE),
    Q3_pAa = quantile(p_Aa, 0.75, na.rm = TRUE),
    aa = median(p_aa, na.rm = TRUE),
    Q1_paa = quantile(p_aa, 0.25, na.rm = TRUE),
    Q3_paa = quantile(p_aa, 0.75, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c("AA", "Aa", "aa"), 
               names_to = "Genotype", 
               values_to = "Median") %>%
  mutate(Q1 = case_when(
    Genotype == "AA" ~ Q1_pAA,
    Genotype == "Aa" ~ Q1_pAa,
    Genotype == "aa" ~ Q1_paa
  ),
  Q3 = case_when(
    Genotype == "AA" ~ Q3_pAA,
    Genotype == "Aa" ~ Q3_pAa,
    Genotype == "aa" ~ Q3_paa
  ))

# Order Factor Levels
df_summary_sep_iqr <- df_summary_sep_iqr %>%
  mutate(Genotype = factor(Genotype, levels = c("AA", "Aa", "aa"))) %>%
  mutate(slope_group = factor(slope_group, levels=c("Positive Slope", "Negative Slope")))

# Step 4: Plot
ggplot(df_summary_sep_iqr, aes(x = Genotype, y = Median, group = cluster, color = factor(cluster))) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +
  labs(x = "Genotype", y = "Value", color = "Cluster") +
  facet_grid(slope_group ~ cluster , scales = "free_y", labeller = labeller(slope_group = c("Positive Slope" = "Positive Slope", "Negative Slope" = "Negative Slope"))) +
  theme_minimal() +
  theme(strip.text.y = element_text(size = 12), 
        strip.text.x = element_text(size = 12), 
        panel.grid.major = element_blank())

# If you wish to strip the background, use: 
# theme(panel.border = element_blank(), 
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank()) 
