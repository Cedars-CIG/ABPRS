library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

##################################################
##
## Part V. BHq selection for both GLM and LM
##
#################################################

######BHq##########
BH<-function(pvalue,q){
  p<-length(pvalue)
  plevel<-(1:p)*q/p
  
  index<-order(pvalue)
  p.seq<-pvalue[index]   
  
  i=1
  index1<-p
  while (i<p) {
    if(p.seq[i]>plevel[i]){
      index1<-i-1
      break
    }
    i<-i+1
  } 
  if(index1>0){
    supp<-index[1:index1]
  }else{
    supp<-NULL
  }
  return(supp) 
}

ORIG_DS_inf_simulation<- function(x, y, family, lasso_est) {
  nn <- length(y)
  pp <- ncol(x)
  x <- cbind(rep(1, nrow(x)), x)
  if(ncol(x) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "gaussian") {
    mu <- as.vector(x%*%lasso_est)
    neg_dloglik_glmnet <- 0 - as.vector(t(x)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- (x*x)/nn
  } else if(family == "binomial") {
    mu <- as.vector(1/(1+exp(-x%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(x)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- (x*(mu*(1-mu)*x))/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(x%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(x)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(x)%*%(mu*x)/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%x
  } else {
    stop("Input family is not supported.")
  }
  
  tau_glmnet <- colSums(neg_ddloglik_glmnet)
  
  theta_glmnet <- 1/tau_glmnet
  
  b_hat_nw <- as.vector(lasso_est - theta_glmnet*neg_dloglik_glmnet)
  se_nw <- theta_glmnet/nn
  pval_nw <- 2*pnorm(abs(b_hat_nw/sqrt(se_nw)), lower.tail=F)
  #se_nw1<-neg_ddloglik_glmnet%*%(neg_ddloglik_glmnet)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw))
}


BH_selection<-function(data_list, family="binomial", alpha=0.1,lambda=3e-3){
  
  prs_train <- data_list$prs_train
  theta_train <- data_list$theta_train
  dat_train <- data_list$dat_train
  
  x<-cbind(prs_train, theta_train)
  
  if(is.null(lambda)){
    fit0<- glm(dat_train ~ x, family = family)
    p_value<-summary(fit0)$coefficients[-c(1,2),4]
  }else{
    fit0 <- glmnet(x, dat_train, family = family,lambda = lambda)
    Theta1 <- as.vector(coef(fit0))
    fit0<-ORIG_DS_inf_simulation(x, y=dat_train, family=family, lasso_est=Theta1)
    p_value<-fit0$pvalue[-c(1,2)]
  }
  selected_snps<-BH(p_value, alpha)
  
  return(list(selected_snps=selected_snps))
}