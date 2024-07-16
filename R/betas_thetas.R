library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

###########################################
##
## Part I. Learning betas and covert theta
##
###########################################

#' Compute Beta Coefficients
#'
#' This function computes beta coefficients and p-values of a generalized 
#' linear model for each SNPs. The model fits the genotype, prs, and covariates 
#' to the phenotype outcome. A binomial distribution should be used for binary
#' outcomes while a gaussian model should be used for continuous outcomes. 
#' 
#' @param pheno a dataframe containing the phenotype (outcome)
#' @param prs a dataframe containing the polygenic risk scores 
#' @param geno a datafrane containing the genotype data whose first column
#' stores the SNP rsID and second column stores the genotype encoding
#' @param cov a dataframe containing the covariates
#' @param family a string naming the regression model used for glm. Use "binomial" 
#' for binary outcomes and "gaussian" for continuous outcomes.
#' @return a dataframe whose columns include the SNP rsID number, beta1, beta2, 
#' and p-value.
#' @export
learning_betas <- function(pheno=pheno_train, geno=snp_train, prs=prs_train, cov=cov_train, family="binomial"){
  p <- ncol(geno)
  betas <- data.frame(SNP=colnames(geno),beta1=rep(NA,p), beta2=rep(NA,p), pval=rep(NA,p))
  geno<-as.matrix(geno)
  for (i in 1:p){
    geno_f <- as.factor(geno[,i])
    if(length(levels(geno_f))>1){
      df <- data.frame(y = pheno, geno_f, prs, cov)
      #Generalized Linear Model
      mod <- glm(y ~ ., data=df, family = family)
      #Retrieve beta1 and beta2
      betas[i,2] <-summary(mod)$coef[2,1]
      betas[i,3] <- summary(mod)$coef[3,1]
      theta <- ifelse(geno[,i]==0 , 0,
                      ifelse(geno[,i] == 1, betas[i,2], betas[i,3]))
      df <- data.frame(y = pheno, theta, prs)
      mod <- glm( y~ ., data=df, family = family)
      betas[i,4] <- summary(mod)$coef[2,4]
      if (i %% 100 == 0) {print(paste("SNP", i, "beta done."))}
    }
  }
  return(betas)
}

#' Convert Theta SNPs
#'
#' This function converts the genotype matrix and the weighted beta coefficients
#' into theta SNPs by replacing homozygous common allele AA with 0, heterozygous 
#' allele Aa with beta1, and homozygous minor allele with beta2.
#'
#' @param geno a dataframe containing genotype data with original SNP encoding
#' @param betas a dataframe containing beta1 and beta2 coefficients (See 
#' \link{learning_betas} to generate the beta dataframe)
#' @return A dataframe with the converted theta SNPs
#' @export
convert_theta <- function(geno, betas){
  theta<-t(t(I(geno==1))*betas$beta1+t(I(geno==2))*betas$beta2)
  theta=as.data.frame(unclass(theta))
  return(theta)
}