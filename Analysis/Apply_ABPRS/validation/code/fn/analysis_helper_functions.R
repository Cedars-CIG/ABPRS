suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cli)
  library(furrr)
})

# m_func -----------------------------------------------------------------------
m_func <- function(
  data,
  .ancestry,
  df_coef,
  rm_cols = c(
    "PMBB_ID", "eMERGE_ID", "ID",
    "IID", "id", "female",
    "ancestry", "y"
    )
) {
  # filter ancestry
  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  data <- data[data$ancestry == .ancestry, ]
  
  y <- data$y
  X <- data |> dplyr::select(-tidyselect::any_of(rm_cols))
  
  df_X_names     <- data.table::data.table(var = names(X))
  df_coef        <- dplyr::inner_join(df_coef, df_X_names, by = "var")
  df_coef_vars   <- df_coef[["var"]]
  df_coef_no_prs <- df_coef_vars[df_coef_vars != "prs"]

  X_all   <- X[, ..df_coef_vars]
  X_prs   <- X_all[, "prs"]
  X_theta <- X_all[, ..df_coef_no_prs]
  
  m1 <- X_prs[[1]] * df_coef[var == "prs", coef]
  m2 <- as.matrix(X_theta) %*% df_coef$coef[-1]
  
  return(list(m1, m2[, 1]))
}

# metric_func ------------------------------------------------------------------
metric_func <- function(
  data,
  .ancestry,
  m1,                # from m_func
  m2,                # from m_func
  split      = TRUE,
  split_seed = 123
) {
  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  data    <- data[data$ancestry == .ancestry, ]
  data$m1 <- m1
  data$m2 <- m2
  data <- data[!is.na(y), ]

  if (nrow(data[, .N, y]) <= 1 || nrow(data) < 100) {

    return(data.table::data.table(ancestry = .ancestry))

  } else {

    if (test_variable_type(data[, y]) == "binary" && sum(data[, y == 1]) < 10) {
      return(data.table::data.table(ancestry = .ancestry))
    }

    if (split == TRUE) {
      data_split <- split_data(data, seed = split_seed)
    } else {
      data_split <- list(
        "train" = data,
        "test"  = data
      )
    }
   
    if (pheno %in% c("AD", "breast_cancer", "hypertension", "t2d")) {
      mod1 <- glm_wc(y ~ prs, data = data_split[["train"]], family = "binomial")
      mod2 <- glm_wc(y ~ m1 + m2, data = data_split[["train"]], family = "binomial")

      mod1_metric_train.ci <- suppressMessages({pROC::ci.auc(data_split[["train"]]$y, predict(mod1$model, newdata = data_split[["train"]], type = "response"))})
      mod2_metric_train.ci <- suppressMessages({pROC::ci.auc(data_split[["train"]]$y, predict(mod2$model, newdata = data_split[["train"]], type = "response"))})

      mod1_metric_train <- mod1_metric_train.ci[[2]]
      mod2_metric_train <- mod2_metric_train.ci[[2]]

      mod1_metric_train_lo <- mod1_metric_train.ci[[1]]
      mod2_metric_train_lo <- mod2_metric_train.ci[[1]]

      mod1_metric_train_hi <- mod1_metric_train.ci[[3]]
      mod2_metric_train_hi <- mod2_metric_train.ci[[3]]

      mod1_metric_train_ci <- print_ci(mod1_metric_train.ci)
      mod2_metric_train_ci <- print_ci(mod2_metric_train.ci)

      mod1_metric.ci <- suppressMessages({pROC::ci.auc(data_split[["test"]]$y, predict(mod1$model, newdata = data_split[["test"]], type = "response"))})
      mod2_metric.ci <- suppressMessages({pROC::ci.auc(data_split[["test"]]$y, predict(mod2$model, newdata = data_split[["test"]], type = "response"))})

      mod1_metric <- mod1_metric.ci[[2]]
      mod2_metric <- mod2_metric.ci[[2]]

      mod1_metric_lo <- mod1_metric.ci[[1]]
      mod2_metric_lo <- mod2_metric.ci[[1]]

      mod1_metric_hi <- mod1_metric.ci[[3]]
      mod2_metric_hi <- mod2_metric.ci[[3]]

      mod1_metric_ci <- print_ci(mod1_metric.ci)
      mod2_metric_ci <- print_ci(mod2_metric.ci)
    } else {
      mod1 <- glm_wc(y ~ prs, data = data_split[["train"]], family = "gaussian")
      mod2 <- glm_wc(y ~ m1 + m2, data = data_split[["train"]], family = "gaussian")
      
      pred1_train       <- predict(mod1$model, newdata = data_split[["train"]])
      mod1_metric_train <- cor(data_split[["train"]]$y, pred1_train)^2
      
      pred2_train       <- predict(mod2$model, newdata = data_split[["train"]])
      mod2_metric_train <- cor(data_split[["train"]]$y, pred2_train)^2

      pred1       <- predict(mod1$model, newdata = data_split[["test"]])
      mod1_metric <- cor(data_split[["test"]]$y, pred1)^2
      mod1_metric_ci <- NA
      pred2       <- predict(mod2$model, newdata = data_split[["test"]])
      mod2_metric <- cor(data_split[["test"]]$y, pred2)^2
      mod2_metric_ci <- NA
    }

    mod1_summary <- summary(mod1$model)$coefficients
    mod2_summary <- summary(mod2$model)$coefficients
    prs_coef   <- mod1_summary["prs", 1]
    prs_coef_p <- mod1_summary["prs", 4]

    m1_coef   <- mod2_summary["m1", 1]
    m1_coef_p <- mod2_summary["m1", 4]
    m2_coef   <- mod2_summary["m2", 1]
    m2_coef_p <- mod2_summary["m2", 4]

     data.table::data.table(
      ancestry                  = .ancestry,
      external_PRS_metric_ci    = mod1_metric_ci,
      AB_PRS_metric_ci          = mod2_metric_ci,
      external_PRS_metric       = mod1_metric,
      AB_PRS_metric             = mod2_metric,
      increase_percent          = (mod2_metric - mod1_metric) / mod1_metric,
      external_PRS_train_metric = mod1_metric_train,
      AB_PRS_train_metric       = mod2_metric_train,
      train_increase_percent    = (mod2_metric_train - mod1_metric_train) / mod1_metric_train,
      prs_coef                  = prs_coef,
      prs_coef_p                = prs_coef_p,
      m1_coef                   = m1_coef,
      m1_coef_p                 = m1_coef_p,
      m2_coef                   = m2_coef,
      m2_coef_p                 = m2_coef_p,
      n_obs_train               = data_split[["train"]][, .N],
      n_obs_test                = data_split[["test"]][, .N],
      test_pct                  = data_split[["test"]][, .N] / data[, .N],
      mod1_warning              = mod1$warning,
      mod2_warning              = mod2$warning
    )
  }
  
}

# metric_rep_func --------------------------------------------------------------
# variation of metric_func that incorporates repetitions of random splits via
# split_seeds argument
metric_rep_func <- function(
  data,
  .ancestry,
  m1,                  # from m_func
  m2,                  # from m_func
  split_seeds = 1:100,
  .cohort
) {
  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  data    <- data[data$ancestry == .ancestry, ]
  data$m1 <- m1
  data$m2 <- m2
  data <- data[!is.na(y), ]

  if (nrow(data[, .N, y]) <= 1 || nrow(data) < 100) {

    return(data.table::data.table(ancestry = .ancestry))

  } else {

    if (test_variable_type(data[, y]) == "binary" && sum(data[, y == 1]) < 10) {
      return(data.table::data.table(ancestry = .ancestry))
    }

    data_splits <- map(split_seeds, \(i) split_data(data, seed = i), .progress = TRUE)

    if (pheno %in% c("AD", "breast_cancer", "hypertension", "t2d")) {
      rep_rslt <- map(
        data_splits,
        \(d) {
          # 1. fit models
            mod1 <- glm_wc(y ~ prs, data = d[["train"]], family = "binomial")
            mod2 <- glm_wc(y ~ m1 + m2, data = d[["train"]], family = "binomial")

          # 2. get predictions
            mod1_pred <- predict(mod1$model, newdata = d[["train"]], type = "response")
            mod2_pred <- predict(mod2$model, newdata = d[["train"]], type = "response")
          
          # 3. get metrics
            mod1_train_metric <- suppressMessages({pROC::auc(d[["train"]]$y, mod1_pred)})[[1]]
            mod2_train_metric <- suppressMessages({pROC::auc(d[["train"]]$y, mod2_pred)})[[1]]

            mod1_metric <- suppressMessages({pROC::auc(d[["test"]]$y, predict(mod1$model, newdata = d[["test"]], type = "response"))})[[1]]
            mod2_metric <- suppressMessages({pROC::auc(d[["test"]]$y, predict(mod2$model, newdata = d[["test"]], type = "response"))})[[1]]

            mod1_summary <- summary(mod1$model)$coefficients
            mod2_summary <- summary(mod2$model)$coefficients
            prs_coef   <- mod1_summary["prs", 1]
            prs_coef_p <- mod1_summary["prs", 4]

            m1_coef   <- mod2_summary["m1", 1]
            m1_coef_p <- mod2_summary["m1", 4]
            m2_coef   <- mod2_summary["m2", 1]
            m2_coef_p <- mod2_summary["m2", 4]

            metrics <- data.table::data.table(
              ancestry                  = .ancestry,
              external_PRS_metric       = mod1_metric,
              AB_PRS_metric             = mod2_metric,
              increase_percent          = (mod2_metric - mod1_metric) / mod1_metric,
              external_PRS_train_metric = mod1_train_metric,
              AB_PRS_train_metric       = mod2_train_metric,
              train_increase_percent    = (mod2_train_metric - mod1_train_metric) / mod1_train_metric
            )
          
          # 4. get stratified ORs
            d_tmp <- data_splits[["train"]] |>
              mutate(pred1 = mod1_pred, pred2 = mod2_pred)

            strat_or(dt = d_tmp, .cohort = cohort)

            
        }, .progress = TRUE
      )
    } else {
      rep_rslt <- map_dfr(
        data_splits,
        \(d) {
          mod1 <- glm_wc(y ~ prs, data = d[["train"]], family = "gaussian")
          mod2 <- glm_wc(y ~ m1 + m2, data = d[["train"]], family = "gaussian")
      
          pred1_train       <- predict(mod1$model, newdata = d[["train"]])
          mod1_train_metric <- cor(d[["train"]]$y, pred1_train)^2
      
          pred2_train       <- predict(mod2$model, newdata = d[["train"]])
          mod2_train_metric <- cor(d[["train"]]$y, pred2_train)^2

          pred1       <- predict(mod1$model, newdata = d[["test"]])
          mod1_metric <- cor(d[["test"]]$y, pred1)^2
          pred2       <- predict(mod2$model, newdata = d[["test"]])
          mod2_metric <- cor(d[["test"]]$y, pred2)^2

          mod1_summary <- summary(mod1$model)$coefficients
          mod2_summary <- summary(mod2$model)$coefficients
          prs_coef   <- mod1_summary["prs", 1]
          prs_coef_p <- mod1_summary["prs", 4]

          m1_coef   <- mod2_summary["m1", 1]
          m1_coef_p <- mod2_summary["m1", 4]
          m2_coef   <- mod2_summary["m2", 1]
          m2_coef_p <- mod2_summary["m2", 4]

          data.table::data.table(
              ancestry                  = .ancestry,
              external_PRS_metric       = mod1_metric,
              AB_PRS_metric             = mod2_metric,
              increase_percent          = (mod2_metric - mod1_metric) / mod1_metric,
              external_PRS_train_metric = mod1_train_metric,
              AB_PRS_train_metric       = mod2_train_metric,
              train_increase_percent    = (mod2_train_metric - mod1_train_metric) / mod1_train_metric
            )
        }, .progress = TRUE
      )
    }
   
    external_PRS_metric <- rep_rslt[, mean(external_PRS_metric)]
    AB_PRS_metric       <- rep_rslt[, mean(AB_PRS_metric)]

    external_PRS_sd <- rep_rslt[, sd(external_PRS_metric)]
    AB_PRS_sd       <- rep_rslt[, sd(AB_PRS_metric)]

    external_PRS_train_metric <- rep_rslt[, mean(external_PRS_train_metric)]
    AB_PRS_train_metric       <- rep_rslt[, mean(AB_PRS_train_metric)]

    external_PRS_train_sd <- rep_rslt[, sd(external_PRS_train_metric)]
    AB_PRS_train_sd       <- rep_rslt[, sd(AB_PRS_train_metric)]

    metric <- data.table(
      ancestry = .ancestry,
      external_PRS_metric = external_PRS_metric,
      external_PRS_sd = external_PRS_sd,
      external_PRS_lo = external_PRS_metric - (qnorm(0.975) * external_PRS_sd / sqrt(length(split_seeds))),
      external_PRS_hi = external_PRS_metric + (qnorm(0.975) * external_PRS_sd/ sqrt(length(split_seeds))), 
      AB_PRS_metric = AB_PRS_metric,
      AB_PRS_sd = AB_PRS_sd,
      AB_PRS_lo = AB_PRS_metric - (qnorm(0.975) * AB_PRS_sd / sqrt(length(split_seeds))),
      AB_PRS_hi = AB_PRS_metric + (qnorm(0.975) * AB_PRS_sd / sqrt(length(split_seeds))),
      increase_percent = (AB_PRS_metric - external_PRS_metric) / external_PRS_metric,
      external_PRS_train_metric = external_PRS_train_metric,
      external_PRS_train_sd = external_PRS_train_sd,
      AB_PRS_train_metric = AB_PRS_train_metric,
      AB_PRS_train_sd = AB_PRS_train_sd,
      increase_percent_train = (AB_PRS_train_metric - external_PRS_train_metric) / external_PRS_train_metric
    )

    strat_or <- map_dfr(
        data_splits,
        \(d) {

        }, .progress = TRUE
    )

  }
  
}

# glm_wc: fit glm model and capture warnings -----------------------------------
glm_wc <- function(
  formula,
  data, 
  family = "binomial"
) {
  warning_occurred <- FALSE  # Initialize warning flag
  
  # Use tryCatch to capture warnings
  model <- tryCatch({
    withCallingHandlers(
      glm(formula = formula, data = data, family = family),
      warning = function(w) {
        warning_occurred <<- TRUE       # Set flag to TRUE if there's a warning
        invokeRestart("muffleWarning")  # Suppress the warning message in the output
      }
    )
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)  # Return NULL if there's an error
  })
  
  # Return both the model and whether a warning occurred
  list(model = model, warning = warning_occurred)
}

# func_SNP_swap ----------------------------------------------------------------
# from Wanling via AB-PRS GitHub
func_SNP_swap <- function(vec) {
  # swap 0 and 2 in .raw
  vec[vec == 0] <- 9
  vec[vec == 2] <- 0
  vec[vec == 9] <- 2
}

# convert_theta ----------------------------------------------------------------
# from Wanling via AB-PRS GitHub
convert_theta <- function(geno, betas_mat) {
  t(t(I(geno == 1)) * betas_mat$beta1 + t(I(geno == 2)) * betas_mat$beta2)
}

# func_process_raw -------------------------------------------------------------
# from Wanling via AB-PRS GitHub
func_process_raw <- function(raw, SNP_list, SNP_complete) {
  ## replace NA with 0
  if (!data.table::is.data.table(raw)) raw <- data.table::as.data.table(raw)
  raw[is.na(raw)] <- 0
  colnames(raw)   <- SNP_complete$SNP
  
  ## check mis-matched effective allele
  SNP_select <- dplyr::inner_join(
      SNP_list,
      SNP_complete,
      by = "SNP"
  )
  
  if (nrow(SNP_select) > 0) {
    SNP_select$match <- (SNP_select$effect_allele == SNP_select$effect_PMBB)
  
    ## select SNPs from .raw
    selected_SNPs <- SNP_select$SNP
    raw_select    <- raw[, ..selected_SNPs]

    for (i in ncol(raw_select)) {
      # if effective allele mis-matched, swap 0 and 2
      if (!SNP_select$match[i]) {
        raw_select[[i]] <- func_SNP_swap(raw_select[[i]])
      }
    }

    if (nrow(SNP_select) == 1) {
      ## if only one SNP selected
      ## force the vector to be column matrix
      raw_select <- as.matrix(raw_select)
      colnames(raw_select) <- SNP_select$SNP[1]
    }

    return(raw_select)
  
  } else {

    return(NULL)

  }
}

# func_prs ---------------------------------------------------------------------
func_prs <- function(rslt) {
  raw  <- rslt[[1]]
  gwas <- rslt[[2]]
  as.matrix(raw) %*% gwas$beta
}

# process_chromosome -----------------------------------------------------------
process_chromosome <- function(i, prs = NULL, GWAS_new, .cohort) {
  # Read the raw file for chromosome i
  raw_dir <- paste0("geno/", .cohort, "/cleaned_chr", i, "_small.raw")
  raw     <- fread(raw_dir)
  
  # Extract IID and subset raw data to keep SNP columns only
  IID <- raw[, 2]          # Save the IID column (if needed for future use)
  raw <- raw[, -c(1:6)] |> # Remove the first 6 columns
    mutate_if(is.numeric, ~round(.))
  
  # Create SNP_complete: SNP names and effect alleles
  SNP_complete <- data.table(
    SNP          = gsub("_.*", "", colnames(raw)),
    effect_cohort = gsub(".*_", "", colnames(raw))
  ) |>
  mutate(
    alt_cohort = substr(SNP, nchar(SNP), nchar(SNP)),
    SNP_long   = paste0(SNP, "_", effect_cohort),
    SNP        = gsub("(:[^:]+){2}$", "", SNP)
  )
  multiallelic <- SNP_complete[, .N, SNP][N > 1, SNP]
  SNP_complete <- SNP_complete[!SNP %in% multiallelic, ]
  these_SNPs   <- SNP_complete[, SNP_long]
  raw          <- raw[, ..these_SNPs]
  SNP_complete <- SNP_complete[, !c("SNP_long")]

  # Process the raw data with the given GWAS
  rslt <- func_process_raw_external_prs(raw, GWAS_new, SNP_complete)
  
  # If there are selected SNPs in chromosome i, update the PRS
  if (nrow(rslt[[2]]) > 0) {
    prs_temp <- func_prs(rslt)
    if (!is.null(prs) && nrow(prs) > 0) {
      prs <- prs + prs_temp
    } else {
      prs <- prs_temp
    }
  }
  return(prs)
}

# func_process_raw_external_prs ------------------------------------------------
func_process_raw_external_prs <- function(raw, .GWAS, SNP_complete) {
  ## replace NA with 0
  raw[is.na(raw)] <- 0
  colnames(raw)   <- SNP_complete$SNP
  
  ## check mis-matched effective allele
  SNP_select       <- inner_join(.GWAS, SNP_complete, by="SNP")
  SNP_select$match <- (SNP_select$effect_allele == SNP_select$effect_cohort)
  
  ## select SNPs from .raw
  selected_SNPs <- SNP_select[, SNP]
  raw_select    <- raw[, ..selected_SNPs]
  
  if (ncol(raw_select) > 0){
    for (i in ncol(raw_select)){
      # if effective allele mis-matched, swap 0 and 2
      if (!SNP_select$match[i]) {
        raw_select[[i]] <- func_SNP_swap(raw_select[[i]])
      }
    }
    return(list(raw_select, SNP_select))
  } else {
    return(list(raw_select, SNP_select))
  }
}

# split_data -------------------------------------------------------------------
split_data <- function(dt, outcome = "y", train_prop = 0.2, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  outcome_type <- test_variable_type(dt[[outcome]])
  
  if (outcome_type == "binary") {
    # Split cases and controls
    dt_cases    <- dt[get(outcome) == 1]
    dt_controls <- dt[get(outcome) == 0]
    
    # Sample 20% of cases and controls for training
    dt_cases[, train := 0]
    dt_cases[sample(.N, .N * train_prop), train := 1]
    
    dt_controls[, train := 0]
    dt_controls[sample(.N, .N * train_prop), train := 1]
    
    # Combine cases and controls back into a single data.table
    dt_split <- rbind(dt_cases, dt_controls)
  } else if (outcome_type == "continuous") {
    dt_split <- dt
    dt_split[, train := 0]
    dt_split[sample(.N, .N * train_prop), train := 1]
  }
  # Create training and testing datasets
  train_data <- dt_split[train == 1]
  test_data  <- dt_split[train == 0]
  
  # Drop the temporary 'train' column
  train_data[, train := NULL]
  test_data[, train := NULL]
  
  # Return list with training and testing data
  list(train = train_data, test = test_data)
}

# split_ids --------------------------------------------------------------------
split_ids <- function(dt, outcome = "y", train_prop = 0.2, seed = 123, id_var = "id") {
  if (!id_var %in% names(dt)) {
    stop(paste0("`id_var` ('", id_var, "') is not in dataset. update `id_var` argument"))
  }
  
  # Set seed for reproducibility
  set.seed(seed)

  outcome_type <- test_variable_type(dt[[outcome]])
  
  if (outcome_type == "binary") {
    # Split cases and controls
    dt_cases    <- dt[get(outcome) == 1]
    dt_controls <- dt[get(outcome) == 0]
    
    # Sample 20% of cases and controls for training
    dt_cases[, train := 0]
    dt_cases[sample(.N, .N * train_prop), train := 1]
    
    dt_controls[, train := 0]
    dt_controls[sample(.N, .N * train_prop), train := 1]
    
    # Combine cases and controls back into a single data.table
    dt_split <- rbind(dt_cases, dt_controls)
  } else if (outcome_type == "continuous") {
    dt_split <- dt
    dt_split[, train := 0]
    dt_split[sample(.N, .N * train_prop), train := 1]
  }
  
  train_ids <- dt_split[train == 1, ][[id_var]]
  test_ids  <- dt_split[train == 0, ][[id_var]]
  
  # Return list with training and testing data
  list(train = train_ids, test = test_ids)
}

# test_variable_type -----------------------------------------------------------
# function to test if a vector is binary or continuous
test_variable_type <- function(vec) {
  # Check if vector is numeric
  if (!is.numeric(vec)) {
    return("non-numeric")
  }
  
  # Remove NAs to avoid interference in checks
  vec_no_na <- na.omit(vec)
  
  # Check for binary indicator (only 0s and 1s)
  if (all(vec_no_na %in% c(0, 1)) && length(unique(vec_no_na)) <= 2) {
    return("binary")
  } else {
    return("continuous")
  }
}

# print_num --------------------------------------------------------------------
print_num <- function(x, r = 3) {
    trimws(format(round(x, r), nsmall = r, big.mark = ","))
}

# print_ci ---------------------------------------------------------------------
print_ci <- function(ci_obj, r = 3, lo_ind = 1, est_ind = 2, hi_ind = 3) {
    paste0(
        print_num(ci_obj[[est_ind]], r = r), " (",
        print_num(ci_obj[[lo_ind]], r = r), ", ",
        print_num(ci_obj[[hi_ind]], r = r), ")"
    )
}

# pred_func --------------------------------------------------------------------
pred_func <- function(data, .ancestry = "EUR", m1, m2) {
  if (!(pheno %in% c("AD", "breast_cancer", "hypertension", "t2d"))) {
    stop("pheno must be one of: AD, breast_cancer, hypertension, t2d")
  }

  # discrete metric: AUC
  # continuous metric: R-squared
  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  data <- data |>
    filter(ancestry == .ancestry) |>
    mutate(m1 = m1, m2 = m2) |>
    filter(!is.na(y))
    data[data$ancestry == .ancestry, ]

  mod1 <- glm_wc(y ~ prs, data = data, family = "binomial")
  mod2 <- glm_wc(y ~ m1 + m2, data = data, family = "binomial")

  data.table(
    id     = data[, ID],
    female = data[, female],
    y      = data[, y],
    pred1  = predict(mod1$model, newdata = data, type = "response"),
    pred2  = predict(mod2$model, newdata = data, type = "response")
  )
}

# OR_func ----------------------------------------------------------------------
OR_func <- function(dt, r1, r0) {
    colnames(dt) <- c("y", "pred_risk")
    binary_output <- classify_by_percentile(dt[["pred_risk"]], range1 = r1, range0 = r0)
    dt <- dt |>
      mutate(X = binary_output)
    
    mod         <- glm(y ~ X, family = "binomial", data = dt)
    odds_ratio  <- exp(coef(mod))[[2]]
    ci          <- suppressMessages({confint(mod)})
    lo          <- exp(ci[2, 1])
    hi          <- exp(ci[2, 2])

    list(
        list(odds_ratio, lo, hi, print_range(r1), print_range(r0), translate_percentile_ranges(list(list(r1, r0)))),
        data.table("or" = odds_ratio, "lo" = lo, "hi" = hi, "range1" = print_range(r1), "range0" = print_range(r0))
    )
}

# strat_or_func ----------------------------------------------------------------
strat_or_func <- function(
  dt,
  pheno_name,
  .cohort = cohort,
  range1 = c(0.9, 1),
  range0 = c(0, 0.9)
) {

    dt_rslt <- data.table(
        sex       = rep(c("Male", "Female", "All"), each = 2),
        model     = rep(c("Pre-trained PRS", "AB-PRS"), 3),
        or        = rep(NA_real_, 6),
        lo        = rep(NA_real_, 6),
        hi        = rep(NA_real_, 6),
        range1    = rep(NA_character_, 6),
        range0    = rep(NA_character_, 6),
        range     = rep(NA_character_, 6)
    )

    # Define a list of conditions and associated rows
    conditions <- list(
        list(filter_sex = 0, pred_col = "pred1", result_row = 1),
        list(filter_sex = 0, pred_col = "pred2", result_row = 2),
        list(filter_sex = 1, pred_col = "pred1", result_row = 3),
        list(filter_sex = 1, pred_col = "pred2", result_row = 4),
        list(filter_sex = NULL, pred_col = "pred1", result_row = 5),
        list(filter_sex = NULL, pred_col = "pred2", result_row = 6)
    )

    # Iterate through conditions
    for (cond in conditions) {
        if (!is.null(cond$filter_sex) && cond$filter_sex == 0 && pheno_name == "breast_cancer")  {
            next 
        } else {
            tmp <- dt |>
                (\(x) { if (!is.null(cond$filter_sex)) filter(x, female == cond$filter_sex) else x})() |>
                select(y, all_of(cond$pred_col))
            tmp_res <- OR_func(tmp, r1 = range1, r0 = range0)[[1]]
            dt_rslt[cond$result_row, 3] <- tmp_res[[1]]
            dt_rslt[cond$result_row, 4] <- tmp_res[[2]]
            dt_rslt[cond$result_row, 5] <- tmp_res[[3]]
            dt_rslt[cond$result_row, 6] <- tmp_res[[4]]
            dt_rslt[cond$result_row, 7] <- tmp_res[[5]]
            dt_rslt[cond$result_row, 8] <- tmp_res[[6]]
        }
    }

    # Exclude Male-specific rows if pheno_name is not "breast_cancer"
    if (pheno_name == "breast_cancer") {
        dt_rslt <- dt_rslt[-c(1, 2), , drop = FALSE]
    }

    return(dt_rslt)
}

# metric_rep_func2 -------------------------------------------------------------
# adaptation of metric_rep_func that incorporates several updates:
#     - allows for arbitrary number of stratified OR estimation within each 
#       split
#     - parallelizes by split to speed up calculation
metric_rep_func2 <- function(
    data,
    .ancestry,
    .id_var = "ID",
    m1,
    m2,
    split_seeds = 1:100,
    .cohort,
    future = TRUE,
    .plan = "multicore",
    n_workers = 16,
    range_list = list(
      list(c(0.9, 1), c(0, 0.9)),
      list(c(0.9, 1), c(0.4, 0.6)),
      list(c(0.9, 1), c(0, 0.1))
    )
) {
   # discrete metric: AUC
  # continuous metric: R-squared
  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  data    <- data[data$ancestry == .ancestry, ]
  data$m1 <- m1
  data$m2 <- m2
  data <- data[!is.na(y), ]

  if (nrow(data[, .N, y]) <= 1 || nrow(data) < 100) {

    return(data.table::data.table(ancestry = .ancestry))

  } else {

    if (test_variable_type(data[, y]) == "binary" && sum(data[, y == 1]) < 10) {
      return(data.table::data.table(ancestry = .ancestry))
    }

    id_splits <- map(split_seeds, \(i) split_ids(data, seed = i, id_var = .id_var), .progress = TRUE)
    # data_splits <- map(split_seeds, \(i) split_data(data, seed = i), .progress = TRUE)

    outcome_type <- test_variable_type(data[["y"]])
    mod_family <- switch(
      outcome_type,
      "binary" = "binomial",
      "continuous" = "gaussian"
    )

    if(future) {
        future::plan(.plan, workers = n_workers)
        tmp_map <- furrr::future_map
        data_obj_size <- round(object.size(data)/1024^2, 2)
        if (data_obj_size > 450) {
          cli_alert_warning("large data object size ({data_obj_size}Mb) requires updating maxSize for parallel operations. setting to {data_obj_size+50}Mb")
          options(future.globals.maxSize = (data_obj_size + 50)*1024^2)
        }
    } else {
        tmp_map <- purrr::map
    }
    
    rep_rslt <- tmp_map(
        id_splits,
        \(d) {
            # 0. split data
            train <- data[get(.id_var) %in% d[["train"]], ]
            test  <- data[get(.id_var) %in% d[["test"]], ]
            
            # 1. fit models
            mod1 <- glm_wc(y ~ prs, data = train, family = mod_family)
            mod2 <- glm_wc(y ~ m1 + m2, data = train, family = mod_family)

            # 2. get predictions
            mod1_pred <- predict(mod1$model, newdata = test, type = "response")
            mod2_pred <- predict(mod2$model, newdata = test, type = "response")
            
            # 3. get metrics
            if (outcome_type == "binary") {
                    mod1_train_metric <- suppressMessages({pROC::auc(train$y, predict(mod1$model, newdata = train, type = "response"))})[[1]]
                    mod2_train_metric <- suppressMessages({pROC::auc(train$y, predict(mod2$model, newdata = train, type = "response"))})[[1]]

                    mod1_metric <- suppressMessages({pROC::auc(test$y, mod1_pred)})[[1]]
                    mod2_metric <- suppressMessages({pROC::auc(test$y, mod2_pred)})[[1]]
            } else if (outcome_type == "continuous") {
                    pred1_train       <- predict(mod1$model, newdata = train)
                    mod1_train_metric <- cor(train$y, pred1_train)^2
                
                    pred2_train       <- predict(mod2$model, newdata = train)
                    mod2_train_metric <- cor(train$y, pred2_train)^2

                    pred1       <- predict(mod1$model, newdata = test)
                    mod1_metric <- cor(test$y, pred1)^2
                    pred2       <- predict(mod2$model, newdata = test)
                    mod2_metric <- cor(test$y, pred2)^2
            }
            
            mod1_summary <- summary(mod1$model)$coefficients
            mod2_summary <- summary(mod2$model)$coefficients
            prs_coef   <- mod1_summary["prs", 1]
            prs_coef_p <- mod1_summary["prs", 4]

            m1_coef   <- mod2_summary["m1", 1]
            m1_coef_p <- mod2_summary["m1", 4]
            m2_coef   <- mod2_summary["m2", 1]
            m2_coef_p <- mod2_summary["m2", 4]

            metrics <- data.table::data.table(
                ancestry                  = .ancestry,
                external_PRS_metric       = mod1_metric,
                AB_PRS_metric             = mod2_metric,
                increase_percent          = (mod2_metric - mod1_metric) / mod1_metric,
                external_PRS_train_metric = mod1_train_metric,
                AB_PRS_train_metric       = mod2_train_metric,
                train_increase_percent    = (mod2_train_metric - mod1_train_metric) / mod1_train_metric
            )
            
            # 4. get stratified ORs (only binary outcomes)
            if (outcome_type == "binary") {
                d_tmp <- test |>
                    mutate(pred1 = mod1_pred, pred2 = mod2_pred)

                strat_or_tabs <- map(
                  range_list,
                  \(r) {
                    strat_or_func(dt = d_tmp, pheno_name = pheno, .cohort = cohort, range1 = r[[1]], range0 = r[[2]])
                  }
                )

                # strat_or_tab <- strat_or_func(dt = d_tmp, pheno_name = pheno, .cohort = cohort)
            } else {
                strat_or_tabs <- list()
            }

            list(
                "metric" = metrics,
                "strat_or" = strat_or_tabs
            )
            
        }, .progress = TRUE
    )

    metrics <- map_dfr(
        seq_along(split_seeds),
        \(i) {
            rep_rslt[[i]][[1]] |> mutate(seed = split_seeds[i])
        }
    )
    strat_ors <- map_dfr(
        seq_along(split_seeds),
        \(i) {
            rbindlist(rep_rslt[[i]][[2]]) |> mutate(seed = split_seeds[i])
        }
    )

    external_PRS_metric <- metrics[, mean(external_PRS_metric)]
    AB_PRS_metric       <- metrics[, mean(AB_PRS_metric)]

    external_PRS_sd <- metrics[, sd(external_PRS_metric)]
    AB_PRS_sd       <- metrics[, sd(AB_PRS_metric)]

    external_PRS_train_metric <- metrics[, mean(external_PRS_train_metric)]
    AB_PRS_train_metric       <- metrics[, mean(AB_PRS_train_metric)]

    external_PRS_train_sd <- metrics[, sd(external_PRS_train_metric)]
    AB_PRS_train_sd       <- metrics[, sd(AB_PRS_train_metric)]

    metric_tab <- data.table(
      ancestry = .ancestry,
      external_PRS_metric = external_PRS_metric,
      external_PRS_sd = external_PRS_sd,
      external_PRS_lo = external_PRS_metric - (qnorm(0.975) * external_PRS_sd / sqrt(length(split_seeds))),
      external_PRS_hi = external_PRS_metric + (qnorm(0.975) * external_PRS_sd/ sqrt(length(split_seeds))), 
      AB_PRS_metric = AB_PRS_metric,
      AB_PRS_sd = AB_PRS_sd,
      AB_PRS_lo = AB_PRS_metric - (qnorm(0.975) * AB_PRS_sd / sqrt(length(split_seeds))),
      AB_PRS_hi = AB_PRS_metric + (qnorm(0.975) * AB_PRS_sd / sqrt(length(split_seeds))),
      increase_percent = (AB_PRS_metric - external_PRS_metric) / external_PRS_metric,
      external_PRS_train_metric = external_PRS_train_metric,
      external_PRS_train_sd = external_PRS_train_sd,
      AB_PRS_train_metric = AB_PRS_train_metric,
      AB_PRS_train_sd = AB_PRS_train_sd,
      increase_percent_train = (AB_PRS_train_metric - external_PRS_train_metric) / external_PRS_train_metric
    )

    return(
        list(
            "metrics_summary" = metric_tab,
            "metrics_tab"     = metrics,
            "strat_or"        = strat_ors
        )
    )
  } 
}

# classify_by_percentile -------------------------------------------------------
classify_by_percentile <- function(
  input_vector,
  range1 = c(0.9, 1),
  range0 = c(0, 0.9)
) {
  if (length(range1) != 2 || length(range0) != 2) {
    stop("Both range1 and range0 should have exactly two elements.")
  }
  
  # Calculate percentile thresholds
  quantiles <- quantile(input_vector, probs = c(range1, range0), na.rm = TRUE)
  
  # Assign binary or NA values
  output <- sapply(input_vector, function(x) {
    if (x >= quantiles[1] && x <= quantiles[2]) {
      return(1)
    } else if (x >= quantiles[3] && x <= quantiles[4]) {
      return(0)
    } else {
      return(NA)
    }
  })
  
  return(output)
}

# print_range ------------------------------------------------------------------
print_range <- function(x) {
  paste0("[", paste0(x, collapse = ", "), "]")
}

# translate_percentile_ranges --------------------------------------------------
translate_percentile_ranges <- function(input_list) {
  sapply(input_list, function(item) {
    top_range   <- item[[1]] # First sub-list represents the 'top'
    other_range <- item[[2]] # Second sub-list represents the 'other'

    # Create top range string
    top_str <- sprintf("top%.1f", diff(top_range))
    
    # Create the other range string
    if (length(other_range) == 2) {
      if (other_range[1] == 0 && other_range[2] < top_range[1]) {
        other_str <- sprintf("bottom%.1f", diff(other_range))
      } else {
        other_str <- sprintf("mid%.1f_%.1f", other_range[1], other_range[2])
      }
    } else {
      stop("Unexpected range format in the second element.")
    }

    # Combine strings
    paste(top_str, other_str, sep = "_")
  })
}
