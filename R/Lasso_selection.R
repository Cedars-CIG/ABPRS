library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

##################################################
##
## Part IV. LASSO selection for both GLM and LM
##
#################################################

#' Lasso Selection
#'
#' This function performs Lasso regression to select important features 
#' (SNPs or others) based on a given dataset and family of models.
#'
#' @param prs_train dataframe containing the training polygenic risk scores
#' @param dat_train dataframe containing the training phenotype data
#' @param theta_train dataframe containing the training thetas
#' @param family specifies the type of model to fit, such as "binomial" for 
#' binary outcomes and "gaussian" for continuous outcomes (default: "binomial")
#' @param lambda regularization parameter for Lasso regression (default: NULL)
#' @return support: a vector indicating indices of selected SNPs from Lasso selection
#' @export
Lasso_selection<-function(prs_train=prs_train, dat_train=dat_train, 
                          theta_train=theta_train, family = "binomial",lambda=NULL){
                          
  if(is.null(lambda)){
    fit0 <- cv.glmnet( as.matrix(cbind(prs_train,theta_train)), dat_train,family = family)
    lambda<-fit0$lambda.min
  }
  fit1 <- glmnet( as.matrix(cbind(prs_train,theta_train)), dat_train, family = family,lambda = lambda)
  Theta <- as.vector(coef(fit1))
  support<-which(Theta[-c(1:2)]!=0)
  
  return(support)
}