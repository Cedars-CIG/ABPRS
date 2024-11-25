# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ROCR)
  library(pROC)
  library(optparse)
  library(furrr)
  library(cli)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

source("fn/analysis_helper_functions.R")
source("fn/liftOver.R")

# Define the options we want to accept
option_list = list(
  make_option(c("-p", "--phenotype"), type = "character", default = NULL, 
              help = "Phenotype outcome", metavar = "character"),
    make_option(c("-c", "--cohort"), type = "character", default = "PMBB", 
              help = "Cohort [default = %default]", metavar = "character")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)
pheno      <- opt$phenotype
cohort     <- tolower(opt$cohort)


# data -------------------------------------------------------------------------
geno_ids <- fread(paste0("data/", tolower(cohort), "_geno_ids"), header = FALSE)[[1]]
prs_dir  <- paste0("external_PRS/", tolower(cohort), "/", pheno, "_external_PRS.txt")
data <- fread(prs_dir) |>
    rename(ID = IID) |>
    filter(ID %in% geno_ids)

for (i in 1:22) {
  theta_dir <- paste0("theta_SNP/", cohort, "/", pheno, "_chr", i, "_theta")
  
  if (file.exists(theta_dir)) {
    theta <- fread(theta_dir)
    data  <- cbind(data, theta)
  }
}

if (tolower(cohort) == "pmbb") {
  pheno_dir <- "data/PMBB_demo.txt"
  demo <- fread(pheno_dir)[, .(PMBB_ID, female)] |> filter(PMBB_ID %in% geno_ids)
  re   <- fread("/static/PMBB/PMBB-Release-2020-2.0/Genotype/PCA/PMBB-Release-2020-2.0_genetic_genotype_ancestries.txt") |>
    select(PMBB_ID, ancestry = Class) |>
    filter(PMBB_ID %in% geno_ids) |>
    mutate(
      ancestry = fifelse(
        ancestry %in% c("UNKNOWN0", "UNKNOWN1", "UNKNOWN2"), "UNKNOWN", ancestry
      )
    )
  demo <- left_join(demo, re, by = "PMBB_ID") |> rename(ID = PMBB_ID)

  if (pheno == "breast_cancer") {
      demo <- demo[female == 1, ]
  }
  phenos    <- fread("data/phecode1.2/PMBB_phenos1.2.txt") |> filter(PMBB_ID %in% geno_ids)
  pheno_dat <- phenos |> select(all_of(c("PMBB_ID", pheno)))
  names(pheno_dat) <- c("ID", "y")

  data <- data[, !duplicated(colnames(data)), with = FALSE]

  demo_pheno <- merge.data.table(demo, pheno_dat, all = FALSE)
  data       <- merge.data.table(demo_pheno, data, all = FALSE)
} else if (tolower(cohort) == "emerge") {
  pheno_dir <- "data/emerge/emerge_demo.txt"
  demo <- fread(pheno_dir)[, .(ID = id, female)]

  ## ancestry
  eur_ids <- fread("/static/eMERGE/eMERGE_III/Info_Files/PCA/kmeans_pca_genetic_ancestry.european", header = FALSE)[[1]]
  afr_ids <- fread("/static/eMERGE/eMERGE_III/Info_Files/PCA/kmeans_pca_genetic_ancestry.african", header = FALSE)[[1]]
  asn_ids <- fread("/static/eMERGE/eMERGE_III/Info_Files/PCA/kmeans_pca_genetic_ancestry.asian", header = FALSE)[[1]]

  eur <- data.table(ID = eur_ids, ancestry = "EUR")
  afr <- data.table(ID = afr_ids, ancestry = "AFR")
  asn <- data.table(ID = asn_ids, ancestry = "ASN")
  re <- bind_rows(eur, afr, asn)

  demo <- left_join(demo, re, by = "ID") |> mutate(ID %in% geno_ids)

  if (pheno == "breast_cancer") {
      demo <- demo[female == 1, ]
  }
  phenos    <- fread("data/emerge/emerge_phenos.txt")
  pheno_dat <- phenos |> select(all_of(c("id", pheno)))
  names(pheno_dat) <- c("ID", "y")

  data <- data[, !duplicated(colnames(data)), with = FALSE]

  demo_pheno <- merge.data.table(demo, pheno_dat)
  data       <- merge.data.table(demo_pheno, data) |> mutate(ID %in% geno_ids)
  data[is.na(ancestry), ancestry := "UNKNOWN"]
} else {
  stop("Cohort must be PMBB or eMERGE")
}


coef_dir <- paste0("UKBB_glm_model/",pheno,"_glm_model.RData")
load(coef_dir)

if (tolower(cohort) == "pmbb") {
  rs_POS_table <- fread("data/rsID_to_chr_pos_table.txt") |>
      rename(var = rsID)
} else if (tolower(cohort) == "emerge") {
  rs_POS_table <- fread("data/rsID_to_chr_pos_table.hg19.txt") |>
    rename(var = rsID)
}


coef_mod2 <- merge.data.table(
    coef_mod2,
    rs_POS_table,
    all.x = TRUE
) |> as.data.table()
coef_mod2[grepl("^rs*", var), var := SNP]
coef_mod2 <- coef_mod2[, .(var, coef)]

rslt    <- m_func(data, "EUR", coef_mod2)
df_rslt <- metric_rep_func2(data, "EUR", m1 = rslt[[1]], m2 = rslt[[2]], .cohort = cohort)

print(df_rslt)

rslt_dir <- paste0("AB_PRS_result/", cohort, "/")
fwrite(df_rslt[["metrics_summary"]], paste0(rslt_dir, pheno, "_V3_result.csv"))
fwrite(df_rslt[["metrics_tab"]], paste0(rslt_dir, pheno, "_V3_rep_result.csv"))
fwrite(df_rslt[["strat_or"]], paste0(rslt_dir, pheno, "_V3_strat_or.csv"))
