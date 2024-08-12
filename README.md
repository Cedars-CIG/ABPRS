# ABPRS

Insert Short Description

## Installation and Loading

1. To install in R library, use:
     ```ruby
     devtools::install_github("oreomilk2005/ABPRS")
     ```
2. To load, use:
     ```ruby
     library(ABPRS)
     ```

## Current Functions:
- ABPRS
- learning_theta_snps
- encoding_theta_snps
- adaptive_variable_selection_binary
- adaptive_variable_selection_continuous
- data_simulation
- model_evaluation


## Basic Function Usage 
     ```ruby
     abprs <- ABPRS(target_prs = target_prs, target_genotype=target_genotype, 
                    training_genotype=training_genotype, validation_genotype=validation_genotype,
                    training_prs=training_prs, validation_prs=validation_prs, 
                    training_phenotype=training_phenotype, validation_phenotype=validation_phenotype,
                    family="binary", biglasso = FALSE, 
                    lam.max = 0.002, lam.min = 6e-05, nlambda = 50,
                    alpha = 0.1, tolerance = 0.025, threshold = 0.01, err = 1e-05,
                    delta = NULL)
     ```

## References

Insert references
