# ABPRS

Insert Short Description

## Installation and Loading

Now: since it's a private package, you will have to download and load the package. 

1. Download the folder on github
2. Then, you can either just load it or install it in your R environment.

     Option1: 
     ```ruby
     devtools::load_all("filepath")
     ```

     Option2:
     ```ruby
     devtools::install("filepath")
     library(ABPRS)
     ```

For Later:

1. To install in R library, use:
     ```ruby
     devtools::install_github("oreomilk2005/ABPRS")
     ```
2. To load, use:
     ```ruby
     library(ABPRS)
     ```

## Current Functions:
- ABPRS: overall function to get weights
- apply_weights: helper function to apply weights
- learning_theta_snps: learn theta1 and theta2 values
- encoding_theta_snps: encode $\theta_{SNPs}$ on genotype
- adaptive_variable_selection: get supports (important $\theta_{SNPs}$)
- data_simulation: generate simulated datasets
- model_evaluation: evaluate different prs models


## Basic Function Usage 
```ruby
weights <- ABPRS(pre_trained_prs = prs_train, validation_prs=prs_val, 
                 training_genotype=geno_train, validation_genotype=geno_val,
                 training_phenotype=pheno_train, validation_phenotype=pheno_val,
                 family="binomial")
```

Please check ABPRS-Tutorial.html for a detailed tutorial. 

## References

Insert references
