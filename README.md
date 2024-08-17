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
- ABPRS
- apply_weights
- learning_theta_snps
- encoding_theta_snps
- adaptive_variable_selection
- data_simulation
- model_evaluation


## Basic Function Usage 
```ruby
weights <- ABPRS(pre_trained_prs = prs_train, validation_prs=prs_val, 
                 training_genotype=geno_train, validation_genotype=geno_val,
                 training_phenotype=pheno_train, validation_phenotype=pheno_val,
                 family="binomial")
```

## References

Insert references
