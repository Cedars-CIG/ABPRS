# AB-PRS

Polygenic risk scores (PRSs) are widely used for predicting genetic risk across
complex diseases and traits, with thousands of pre-trained models constructed
using diverse methodologies and data sources. Despite their widespread utility,
few approaches leverage these pre-trained PRS to further improve their predictive
performance. We present Adaptively Boosting of Pre-trained PRS (AB-PRS), a
fine-tuning framework that refines pre-trained PRS as baseline model through
adaptive variable selection and model boosting to uncover signals missed by the
original models. Simulations show AB-PRS effectively identifies signals orthogo-
nal to pre-trained PRS while controlling false discovery rates. Using UK Biobank
data, we fine-tuned pre-trained PRS for binary diseases and continuous traits, val-
idating results across three independent datasets: All of Us, eMERGE, and Penn
Medicine Biobank. Real data analyses confirmed improved model prediction, risk
stratification, and robustness. AB-PRS provides a scalable solution to enhance
existing and future PRS, advancing predictive capabilities in diverse datasets.

## Installation and Loading

1. To install in R library, use:
     ```ruby
     devtools::install_github("Cedars-CIG/ABPRS")
     ```
2. To load, use:
     ```ruby
     library(ABPRS)
     ```

## Tutorial
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

To Be Added
