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
     devtools::install_github("wuoli/ABPRS")
     ```
2. To load, use:
     ```ruby
     library(ABPRS)
     ```
     
# Tutorial
Please refer to https://cedars-cig.github.io/ABPRS/ for a getting-started tutorial. 

## Current Functions:
- ABPRS: overall function to get weights
- apply_weights: helper function to apply weights
- learning_theta_snps: learn theta1 and theta2 values
- encoding_theta_snps: encode $\theta_{SNPs}$ on genotype
- adaptive_variable_selection: get supports (important $\theta_{SNPs}$)
- data_simulation: generate simulated datasets
- model_evaluation: evaluate various PRS models


## Basic Function Usage 
```ruby
weights <- ABPRS(pre_trained_prs = training_prs, validation_prs=validation_prs, 
                 training_genotype=training_genotype, validation_genotype=validation_genotype,
                 training_phenotype=training_phenotype, validation_phenotype=validation_phenotype,
                 family="binomial", covariate=NULL, biglasso=FALSE, 
                 lam.max=NULL, lam.min=NULL, nlambda=100,
                 alpha=0.1, tolerance=0.025, threshold=0.01,
                 err=1e-5, delta=NULL)
```

## References

Insert references
