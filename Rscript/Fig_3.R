library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(latex2exp)

cohort <- "UKBB" # option: UKBB, emerge, aou, pmbb
args <- c(paste0("metric_",cohort,".csv"), paste0("metric_plot_",cohort,".pdf"))

df_rslt <- read.csv(args[1], header=T)

## ---- read data ----
## Prefer the uploaded CSV if present; otherwise fall back to metric_UKBB.csv

df_rslt <- read.csv(args[1], header = TRUE)

## ensure rep exists for paired test
if (!("rep" %in% names(df_rslt))) {
  df_rslt$rep <- ave(df_rslt$pheno, df_rslt$source, FUN = seq_along)
}

########################################
## helper: mean + CI for bars
########################################
func <- function(vec, CI_level = 0.95){
  m <- mean(vec)
  s <- sd(vec)
  q <- qnorm(1 - 0.5 * (1 - CI_level))
  lw <- m - q * s / sqrt(length(vec))
  up <- m + q * s / sqrt(length(vec))
  list(m, lw, up)
}

########################################
## build summary table used for bars
########################################
make_summary <- function(df_rslt, source_name, dataset_label, phenotype_vec, cont_disc){
  df1 <- data.frame(
    Phenotype = phenotype_vec,
    Model = rep("Pre-trained PRS", length(phenotype_vec)),
    Metric = rep(NA_real_, length(phenotype_vec)),
    lw = rep(NA_real_, length(phenotype_vec)),
    up = rep(NA_real_, length(phenotype_vec))
  )
  df2 <- data.frame(
    Phenotype = phenotype_vec,
    Model = rep("AB-PRS", length(phenotype_vec)),
    Metric = rep(NA_real_, length(phenotype_vec)),
    lw = rep(NA_real_, length(phenotype_vec)),
    up = rep(NA_real_, length(phenotype_vec))
  )
  
  for (i in seq_along(phenotype_vec)) {
    ph <- phenotype_vec[i]
    dfi <- df_rslt[df_rslt$pheno == ph & df_rslt$source == source_name, ]
    if (nrow(dfi) == 0) next
    
    rslt1 <- func(dfi$PRS)
    df1$Metric[i] <- rslt1[[1]]
    df1$lw[i] <- rslt1[[2]]
    df1$up[i] <- rslt1[[3]]
    
    rslt2 <- func(dfi$AB_PRS)
    df2$Metric[i] <- rslt2[[1]]
    df2$lw[i] <- rslt2[[2]]
    df2$up[i] <- rslt2[[3]]
  }
  
  out <- rbind(df1, df2)
  out$dataset <- dataset_label
  out$pheno <- cont_disc
  out
}

## consortium / pgs have 8 phenotypes; finngen has 5 (only bmi for continuous)
phenotype_vec_8 <- c("hypertension","t2d","AD","breast_cancer","ldl","hdl","bmi","cholesterol")
cont_disc_8 <- c("Discrete Phenotype","Discrete Phenotype","Discrete Phenotype","Discrete Phenotype",
                 "Continuous Phenotype","Continuous Phenotype","Continuous Phenotype","Continuous Phenotype")

phenotype_vec_fg <- c("hypertension","t2d","AD","breast_cancer","bmi")
cont_disc_fg <- c("Discrete Phenotype","Discrete Phenotype","Discrete Phenotype","Discrete Phenotype",
                  "Continuous Phenotype")

df_long_consortium <- make_summary(df_rslt, "consortium", "Consortium", phenotype_vec_8, cont_disc_8)
df_long_finngen    <- make_summary(df_rslt, "finngen", "FinnGen", phenotype_vec_fg, cont_disc_fg)
df_long_pgs        <- make_summary(df_rslt, "pgs", "PGS-Catalog", phenotype_vec_8, cont_disc_8)

df_long_all <- rbind(df_long_consortium, df_long_finngen, df_long_pgs)

########################################
## p-values: per (dataset, pheno_type, Phenotype)
########################################
pheno_type_map <- c(
  hypertension  = "Discrete Phenotype",
  t2d           = "Discrete Phenotype",
  AD            = "Discrete Phenotype",
  breast_cancer = "Discrete Phenotype",
  ldl           = "Continuous Phenotype",
  hdl           = "Continuous Phenotype",
  bmi           = "Continuous Phenotype",
  cholesterol   = "Continuous Phenotype"
)

########################################
## raw data points used for dot overlays
########################################
df_point_all <- df_rslt %>%
  mutate(
    dataset = case_when(
      source == "pgs"        ~ "PGS-Catalog",
      source == "consortium" ~ "Consortium",
      source == "finngen"    ~ "FinnGen",
      TRUE ~ as.character(source)
    ),
    pheno_type = unname(pheno_type_map[pheno])
  ) %>%
  filter(!is.na(pheno_type)) %>%
  pivot_longer(
    cols = c(PRS, AB_PRS),
    names_to = "Model_raw",
    values_to = "Metric"
  ) %>%
  mutate(
    Phenotype = pheno,
    Model = case_when(
      Model_raw == "PRS"    ~ "Pre-trained PRS",
      Model_raw == "AB_PRS" ~ "AB-PRS",
      TRUE ~ Model_raw
    ),
    pheno = pheno_type
  ) %>%
  select(Phenotype, Model, Metric, dataset, pheno, rep)

df_p <- df_rslt %>%
  mutate(
    dataset = case_when(
      source == "pgs"        ~ "PGS-Catalog",
      source == "consortium" ~ "Consortium",
      source == "finngen"    ~ "FinnGen",
      TRUE ~ as.character(source)
    ),
    pheno_type = unname(pheno_type_map[pheno])
  ) %>%
  filter(!is.na(pheno_type)) %>%
  group_by(dataset, pheno_type, pheno) %>%
  summarise(
    p = tryCatch(t.test(PRS, AB_PRS, paired = TRUE)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      is.na(p)   ~ "\u2248",
      p < 0.001  ~ "***",
      p < 0.01   ~ "**",
      p < 0.05   ~ "*",
      TRUE       ~ "\u2248"
    ),
    p_label = case_when(
      is.na(p)  ~ NA_character_,
      p < 0.001 ~ paste0("p = ", formatC(p, format = "e", digits = 3)),
      TRUE      ~ paste0("p = ", formatC(p, format = "f", digits = 3))
    )
  )

## y-position based on bar upper CI (so brackets sit above error bars)
y_base <- df_long_all %>%
  group_by(dataset, pheno = pheno, Phenotype) %>%
  summarise(y_base = max(up, na.rm = TRUE), .groups = "drop")

df_p <- df_p %>%
  left_join(
    y_base %>% rename(pheno_type = pheno),
    by = c("dataset", "pheno_type", "pheno" = "Phenotype")
  ) %>%
  mutate(
    y = y_base + 0.01,          # bump; increase if needed
    y.tick = y - 0.003          # tick length
  )

########################################
## Plot settings (controlled dodge for bracket x-positions)
########################################
DODGE_W <- 0.8
OFFSET <- DODGE_W / 4  # for 2 groups, centers at x +/- w/4

data_labs <- c("Consortium" = "Consortium", "FinnGen" = "FinnGen", "PGS-Catalog" = "PGS-Catalog")

## Helper to build bracket data for a given phenotype ordering + panel subset
make_brackets <- function(df_p_sub, phen_levels){
  if (nrow(df_p_sub) == 0) return(df_p_sub)
  
  df_p_sub <- df_p_sub %>%
    mutate(
      Phenotype = factor(pheno, levels = phen_levels),
      x = as.numeric(Phenotype),
      x1 = x - OFFSET,
      x2 = x + OFFSET
    )
  
  df_p_sub
}

########################################
## ---- p1: Discrete Phenotype ----
########################################
disc_levels <- c("AD", "breast_cancer", "hypertension", "t2d")

df_disc <- df_long_all %>%
  filter(pheno == "Discrete Phenotype") %>%
  mutate(
    dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
    Model = factor(Model, levels = c("Pre-trained PRS", "AB-PRS")),
    Phenotype = factor(Phenotype, levels = disc_levels)
  )

df_disc_point <- df_point_all %>%
  filter(pheno == "Discrete Phenotype") %>%
  mutate(
    dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
    Model = factor(Model, levels = c("Pre-trained PRS", "AB-PRS")),
    Phenotype = factor(Phenotype, levels = disc_levels)
  )

p1 <- ggplot(df_disc, aes(x = Phenotype, y = Metric, fill = Model)) +
  geom_col(position = position_dodge(width = DODGE_W), width = 0.7, color = NA) +
  geom_point(
    data = df_disc_point,
    aes(x = Phenotype, y = Metric, group = Model),
    position = position_dodge(width = DODGE_W),
    shape = 16, size = 1.7, color = "black", alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = lw, ymax = up),
    position = position_dodge(width = DODGE_W),
    width = 0.18, linewidth = 0.8, color = "black"
  ) +
  scale_fill_manual(values = c("#fccccb", "#b0d992")) +
  coord_cartesian(ylim = c(0.5, NA)) +
  labs(y = "AUC") +
  facet_grid(pheno ~ dataset, labeller = labeller(dataset = data_labs)) +
  scale_x_discrete(labels = c("AD" = "AD", "hypertension" = "Hypertension",
                              "breast_cancer" = "Breast\nCancer", "t2d" = "T2D")) +
  theme(
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    text = element_text(size = 16),
    strip.text = element_text(size = 17),
    axis.text.x = element_text(size = 15)
  )

## brackets for discrete
br_disc <- df_p %>%
  filter(pheno_type == "Discrete Phenotype") %>%
  make_brackets(disc_levels) %>%
  mutate(dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
         pheno = "Discrete Phenotype")

p1 <- p1 +
  geom_segment(data = br_disc, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = br_disc, aes(x = x1, xend = x1, y = y.tick, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = br_disc, aes(x = x2, xend = x2, y = y.tick, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = br_disc, aes(x = x, y = y + 0.005, label = label),
            inherit.aes = FALSE, size = 5) +
  geom_text(data = br_disc, aes(x = x, y = y + 0.018, label = p_label),
            inherit.aes = FALSE, size = 3.5)

########################################
## ---- p2: Continuous Phenotype ----
########################################
cont_levels <- c("bmi", "cholesterol", "hdl", "ldl")

df_cont <- df_long_all %>%
  filter(pheno == "Continuous Phenotype") %>%
  mutate(
    dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
    Model = factor(Model, levels = c("Pre-trained PRS", "AB-PRS")),
    Phenotype = factor(Phenotype, levels = cont_levels)
  )

df_cont_point <- df_point_all %>%
  filter(pheno == "Continuous Phenotype") %>%
  mutate(
    dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
    Model = factor(Model, levels = c("Pre-trained PRS", "AB-PRS")),
    Phenotype = factor(Phenotype, levels = cont_levels)
  )

p2 <- ggplot(df_cont, aes(x = Phenotype, y = Metric, fill = Model)) +
  geom_col(position = position_dodge(width = DODGE_W), width = 0.7, color = NA) +
  geom_point(
    data = df_cont_point,
    aes(x = Phenotype, y = Metric, group = Model),
    position = position_dodge(width = DODGE_W),
    shape = 16, size = 1.7, color = "black", alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = lw, ymax = up),
    position = position_dodge(width = DODGE_W),
    width = 0.18, linewidth = 0.8, color = "black"
  ) +
  scale_fill_manual(values = c("#fccccb", "#b0d992")) +
  labs(y = TeX("$R^2$")) +
  facet_grid(pheno ~ dataset, labeller = labeller(dataset = data_labs)) +
  scale_x_discrete(labels = c("bmi" = "BMI", "cholesterol" = "Cholesterol",
                              "hdl" = "HDL", "ldl" = "LDL")) +
  theme(
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 16),
    strip.text = element_text(size = 17),
    axis.text.x = element_text(size = 15)
  )

## brackets for continuous
br_cont <- df_p %>%
  filter(pheno_type == "Continuous Phenotype") %>%
  make_brackets(cont_levels) %>%
  mutate(dataset = factor(dataset, levels = c("PGS-Catalog", "Consortium", "FinnGen")),
         pheno = "Continuous Phenotype")

p2 <- p2 +
  geom_segment(data = br_cont, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = br_cont, aes(x = x1, xend = x1, y = y.tick, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = br_cont, aes(x = x2, xend = x2, y = y.tick, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = br_cont, aes(x = x, y = y + 0.005, label = label),
            inherit.aes = FALSE, size = 5) +
  geom_text(data = br_cont, aes(x = x, y = y + 0.018, label = p_label),
            inherit.aes = FALSE, size = 3.5)

########################################
## combine
########################################
p <- p1 + p2 + plot_layout(nrow = 2)

## optional save
ggsave(args[2], p, device = cairo_pdf, width = 45, height = 20, units = "cm")