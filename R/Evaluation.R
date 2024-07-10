library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

###################################################
##
## Part II. AUC test for GLM and predict loss for LM 
##
###################################################

#Expit function (Logistic Sigmoid Function)
expit<-function(x){
  1/(1+exp(-x))
}

predict_test<-function(x,est){
  mu<-est[1]+x%*%est[-1]
  return(expit(mu))
}

#Refined Debiased Lasso
REF_DS_inf <- function(x, y, family, lasso_est, delta) {
  nn <- length(y)
  X <- cbind(1, x)
  p<-ncol(X)
  if(p != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.") 
  }
  if(family == "gaussian") {
    mu <- as.vector(X%*%lasso_est)
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- (t(X)%*%X)/nn
  }else if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%((mu*(1-mu))*X)/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
  } else {
    stop("Input family is not supported.")
  }
  
  theta_inv <-tryCatch({
    solve(neg_ddloglik_glmnet)
  },error=function(e){ tryCatch({
    solve(neg_ddloglik_glmnet+delta*diag(rep(1,p)))
  },error=function(e){
    cat("delta is too small \n")
  })
  }
  )
  b_hat_inv <- as.vector(lasso_est - theta_inv%*%neg_dloglik_glmnet)
  se_inv <- sqrt(diag(theta_inv))/sqrt(nn)
  pval_inv <- 2*pnorm(abs(b_hat_inv/se_inv), lower.tail=F)
  
  return(list(est=b_hat_inv, se=se_inv, pvalue=pval_inv, theta=theta_inv))
}

#' Area Under the Curve Score for Binary Outcomes
#'
#' This function computes the Area Under the Curve (AUC) scores for the 
#' test data.
#' 
#' @param support vector of indices specifying which columns of thetas to
#' include in the model
#' @param prs_train dataframe containing the training polygenic risk scores
#' @param prs_test dataframe containing the testing polygenic risk scores
#' @param dat_train dataframe containing the training phenotype data
#' @param dat_test dataframe containing the testing phenotype data
#' @param theta_train dataframe containing the training thetas
#' @param theta_test dataframe containing the testing thetas
#' @param lasso boolean flag to determine whether to use lasso regression
#' @param lambda tuning parameters for lasso regression
#' @param delta additional tuning parameter for the Lasso regression; used for
#' refined debiased lasso regression (default: 0)
#' @return the AUC score of the data's prediction performance
#' @export
AUC_test<-function(support, prs_train=prs_train, prs_test=prs_test, dat_train=dat_train,
                   dat_test=dat_test, theta_train=theta_train, theta_test=theta_test,
                   lasso=FALSE, lambda=NULL,delta=0){
  
  #If there is support, retrieve support from thetas; if not, use all data
  if(length(support)>0){
    df_train_2 <- data.frame(prs=prs_train, theta_new = theta_train[,support])
    df_test_2 <- data.frame(prs=prs_test, theta_new = theta_test[,support])
  }else{
    df_train_2 <- data.frame(prs=prs_train)
    df_test_2 <- data.frame(prs=prs_test)
  }
  
  # If lasso, fit data using lasso regression model and produce performance metric;
  # if not lasso, directly calculate performance metric using binomial regression
  if(lasso){
    df_train_2<-as.matrix(df_train_2)
    if(is.null(lambda)){
      fit0 <- cv.glmnet(df_train_2, dat_train,nfolds = 5, family = "binomial")
      lambda<-fit0$lambda.min
    }
    fit0 <- glmnet(df_train_2, dat_train, family = "binomial",lambda = lambda)
    Theta <- as.vector(coef(fit0))
    fit0<-REF_DS_inf(x=df_train_2, y=dat_train, family="binomial", lasso_est=Theta, delta)
    Theta<-as.vector(fit0$est)
    
    prediction <- prediction(predict_test(as.matrix(df_test_2),Theta), dat_test)
    mod2_val_auc <- performance(prediction, measure="auc")@y.values[[1]]
  }else{
    # Binomial regression
    mod2 <- glm(dat_train ~ ., data = df_train_2, family="binomial")
    # Obtain predictions
    prediction <- prediction(predict(mod2, newdata=df_test_2, type="response"), dat_test)
    # Find auc score from prediction
    mod2_val_auc <- performance(prediction, measure="auc")@y.values[[1]]
  }
  return(mod2_val_auc)
  
}

#' Predict loss for LM on test data
#'
#' This function computes the mean squared error between the predicted and
#' actual values.
#' 
#' @param support vector of indices specifying which columns of thetas to
#' include in the model
#' @param prs_train dataframe containing the training polygenic risk scores
#' @param prs_test dataframe containing the testing polygenic risk scores
#' @param dat_train dataframe containing the training phenotype data
#' @param dat_test dataframe containing the testing phenotype data
#' @param theta_train dataframe containing the training thetas
#' @param theta_test dataframe containing the testing thetas
#' @param lasso boolean flag to determine whether to use lasso regression
#' @param lambda tuning parameters for lasso regression (default: NULL)
#' @param delta additional tuning parameter for the Lasso regression; used for
#' refined debiased lasso regression (default: 0)
#' @return Mean squared error (MSE) score 
#' @export
Loss_test<-function(support, prs_train=prs_train, prs_test=prs_test, dat_train=dat_train,
                    dat_test=dat_test, theta_train=theta_train, theta_test=theta_test,
                    lasso=FALSE, lambda=NULL, delta=0){

  if(length(support)>0){
    df_train_2 <- data.frame(prs=prs_train, theta_new = theta_train[,support])
    df_test_2 <- data.frame(prs=prs_test, theta_new = theta_test[,support])
  }else{
    df_train_2 <- data.frame(prs=prs_train)
    df_test_2 <- data.frame(prs=prs_test)
  }
  
  
  if(lasso){
    df_train_2<-as.matrix(df_train_2)
    if(is.null(lambda)){
      fit0 <- cv.glmnet(df_train_2, dat_train,nfolds = 5, family = "gaussian")
      lambda<-fit0$lambda.min
    }
    fit0 <- glmnet(df_train_2, dat_train, family = "gaussian",lambda = lambda)
    Theta <- as.vector(coef(fit0))
    fit0<-REF_DS_inf(x=df_train_2, y=dat_train, family="gaussian", lasso_est=Theta, delta)
    Theta<-as.vector(fit0$est)
    
    prediction <- (as.matrix(cbind(1,df_test_2))%*%Theta-dat_test)^2
    mod2_val_auc <- mean(prediction)
    
  }else{
    
    mod2 <- glm(dat_train ~ ., data = df_train_2, family="gaussian")
    # training AUC diff
    prediction <- (predict(mod2, newdata=df_test_2, type="response")- dat_test)^2
    mod2_val_auc <- mean(prediction)
  }
  return(mod2_val_auc)
}