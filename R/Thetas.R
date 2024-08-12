library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

#' Compute Theta SNPs for ABPRS
#'
#' This function computes theta SNPs and their corresponding p-values for each SNP.
#' Specifically, the function fits a regression model to the phenotype outcome using the
#' genotype, PRS, and covariates. The function uses a binomial model for binary
#' outcomes and a Gaussian model for continuous outcomes.
#' 
#' @param phenotype a dataframe containing the phenotype
#' @param prs a dataframe containing the polygenic risk scores 
#' @param genotype a dataframe containing the genotype matrix with additive encoding
#' whose columns represent SNPs and rows represent individuals.
#' @param covariate a dataframe containing covariates (default: NULL)
#' @param family a string describing the category of phenotype in the dataset. 
#' Use "binary" for binary outcomes and "continuous" for continuous outcomes
#' (default: "binary").
#' @return A dataframe whose columns include the SNP rsID, theta1 (theta encoding 
#' for heterozygous alleles Aa), theta2 (theta encoding for homozygous alternative 
#' alleles aa), and p-value.
#' @export
learning_theta_snps <- function(phenotype, genotype, prs, covariate=NULL, family="binary"){
  
  if(family == "binary"){
    family <- "binomial"
  }else if(family == "continuous"){
    family <- "gaussian"
  }
  
  p <- ncol(genotype)
  if(is.null(covariate)){
    covariate <- matrix(nrow=length(phenotype), ncol=1, 0)
  }
  thetas <- data.frame(SNP=colnames(genotype),theta1=rep(NA,p), theta2=rep(NA,p), pval=rep(NA,p))
  genotype<-as.matrix(genotype)
  for (i in 1:p){
    geno_f <- as.factor(genotype[,i])
    if(length(levels(geno_f))>1){
      df <- data.frame(y = phenotype, geno_f, prs, covariate)
      #Generalized Linear Model
      mod <- glm(y ~ ., data=df, family = family)
      #Retrieve theta1 and theta2
      thetas[i,2] <-summary(mod)$coef[2,1]
      thetas[i,3] <- summary(mod)$coef[3,1]
      theta <- ifelse(genotype[,i]==0 , 0,
                      ifelse(genotype[,i] == 1, thetas[i,2], thetas[i,3]))
      df <- data.frame(y = phenotype, theta, prs)
      mod <- glm( y~ ., data=df, family = family)
      thetas[i,4] <- summary(mod)$coef[2,4]
      #if (i %% 100 == 0) {print(paste("SNP", i, "theta done."))}
    }
  }
  return(thetas)
}

#' Convert Additive Encoding to Theta SNP Encoding
#'
#' This function converts a genotype matrix with additive encoding into a
#' theta-encoded matrix. It replaces the encoding of heterozygous alleles (Aa) 
#' with theta1 and homozygous alternative alleles (aa) with theta2, while retaining 
#' the encoding of homozygous reference alleles (AA) as 0.
#'
#' @param genotype a dataframe containing genotype data with additive SNP encoding.
#' @param thetas a dataframe containing theta encoding for genotype Aa (theta1) 
#' and aa (theta2) (See \link{learning_theta_snps} for the function to generate 
#' this dataframe).
#' @return A dataframe containing the SNPs with theta encoding. 
#' @export
encoding_theta_snps <- function(genotype, thetas){
  encoded_theta<-t(t(I(genotype==1))*thetas$theta1+t(I(genotype==2))*thetas$theta2)
  encoded_theta<-as.data.frame(unclass(encoded_theta))
  return(encoded_theta)
}