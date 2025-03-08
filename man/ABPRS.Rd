% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABPRS.R
\name{ABPRS}
\alias{ABPRS}
\title{ABPRS Function}
\usage{
ABPRS(
  pre_trained_prs,
  validation_prs,
  training_phenotype,
  validation_phenotype,
  training_genotype,
  validation_genotype,
  family,
  covariate = NULL,
  biglasso = FALSE,
  lam.max = NULL,
  lam.min = NULL,
  nlambda = 100,
  alpha = 0.1,
  tolerance = 0.025,
  threshold = 0.01,
  err = 1e-05,
  delta = NULL
)
}
\arguments{
\item{pre_trained_prs}{dataframe containing the pre-trained training polygenic risk scores}

\item{validation_prs}{dataframe containing the pre-trained validation polygenic risk scores}

\item{training_phenotype}{vector containing the training phenotype data}

\item{validation_phenotype}{vector containing the validation phenotype data}

\item{training_genotype}{dataframe containing the training genotype data}

\item{validation_genotype}{dataframe containing the validation genotype data}

\item{family}{string specifying the type of model to fit a \link[stats]{glm}.
Use "binomial" for binary outcomes and "gaussian" for continuous outcomes.}

\item{covariate}{dataframe containing the covariates (default: NULL)}

\item{biglasso}{logical flag indicating whether to use big data optimization.
If true, then the program uses \link[biglasso]{biglasso}. Or else, it uses
the standard \link[glmnet]{glmnet} (default: FALSE)}

\item{lam.max}{maximum value of the regularization parameter (lambda) to be
considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. When NULL (default),
glmnet automatically generates the lambda sequence.}

\item{lam.min}{minimum value of the regularization parameter (lambda) to be
considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. When NULL (default),
glmnet automatically generates the lambda sequence.}

\item{nlambda}{number of different lambda values to be evaluated between
lam.max and lam.min. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}.}

\item{alpha}{desired FDR control level.}

\item{tolerance}{tolerance level for noise in mirror statistics.}

\item{threshold}{threshold value ``v” for determining the significance of
difference between training improvement and validation improvement,
see equation (3) of ABPRS paper.}

\item{err}{error threshold for model performance improvement.}

\item{delta}{if ``NULL” then using GLM to fit the regression model on the validation dataset (default).
If a small value is specified, the debiased lasso method is applied, where the Hessian matrix is adjusted
by adding an identity matrix weighted by delta to ensure invertibility.}
}
\value{
A dataframe with three columns. The first column named "Weight" stores
the weights for pre-trained PRS scores and selected \eqn{\theta_{SNPs}}. The
second and third columns stores the theta1 (theta encoding for heterozygous
alleles Aa) and theta2 (theta encoding for homozygous alternative alleles aa)
of the selected \eqn{\theta_{SNPs}}, which are useful for encoding the genotype
when producing adaptively-boosted polygenic risk scores (AB-PRSs).
}
\description{
This function runs the ABPRS method and returns the weights of
the pre-trained prs and selected \eqn{\theta_{SNPs}} used for ABPRS calculation.
For the selected \eqn{\theta_{SNPs}}, it also returns their corresponding
theta1 (theta encoding for heterozygous alleles Aa) and theta2 (theta encoding
for homozygous alternative alleles aa) values, which are important for encoding
the genotype matrices when producing adaptively-boosted polygenic risk scores (AB-PRSs).
}
