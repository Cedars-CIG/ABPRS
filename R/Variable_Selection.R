library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

# Helper Functions
FDR_control<-function(beta1, beta2, q){
  #Mirror Statistics
  M<- sign(beta1*beta2)*(abs(beta1)*I(abs(beta1)<abs(beta2))+abs(beta2)*I(abs(beta1)>=abs(beta2)))
  tau.seq<-abs(M[I(M!=0)])
  #Determine threshold value tau
  tau.seq<-tau.seq[order(tau.seq)]
  for (tau in c(0,tau.seq)) {
    fdp<- sum(I(M<=(-tau)))/max(sum(I(M>=tau)),1)
    if(fdp<=q){
      tau.hat<-tau
      break
    } 
    tau.hat<-tau
  }
  #Get support
  p<-length(M)
  index<-which(M>=tau.hat)
  S<-rep(0,p)
  S[index]<-1
  return(S)
}

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

#' Adaptive Variable Selection for Binary Outcomes
#'
#' This function identifies important \eqn{\theta_{SNPs}} through adaptive variable
#' selection and FDR control for binary outcomes. 
#'
#' @param training_prs dataframe containing the training polygenic risk scores
#' @param validation_prs dataframe containing the validation polygenic risk scores
#' @param training_phenotype dataframe containing the training phenotype data
#' @param validation_phenotype dataframe containing the validation phenotype data
#' @param training_theta_encoding dataframe containing the training thetas
#' @param validation_theta_encoding dataframe containing the validation thetas
#' @param biglasso logical flag indicating whether to use big data optimization. 
#' If true, then the program uses \link[biglasso]{biglasso}. Or else, it uses 
#' the standard \link[glmnet]{glmnet} (default: FALSE)
#' @param lam.max maximum value of the regularization parameter (lambda) to be 
#' considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param lam.min minimum value of the regularization parameter (lambda) to be 
#' considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param nlambda number of different lambda values to be evaluated between 
#' lam.max and lam.min. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param alpha desired FDR control level.
#' @param tolerance tolerance level for noise in mirror statistics.
#' @param threshold threshold value for determining significance in mirror 
#' statistics.
#' @param err error threshold for model performance improvement.
#' @param delta additional tuning parameter for the FDR control.
#' @return A vector of selected \eqn{\theta_{SNPs}}. 
#' @export
adaptive_variable_selection_binary<-function(training_prs, validation_prs, 
                                             training_phenotype, validation_phenotype, 
                                             training_theta_encoding, validation_theta_encoding,
                                             biglasso=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                                             alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-5, delta=NULL){
  
  p<-ncol(training_theta_encoding)
  
  # for training
  df_train_1 <- data.frame(prs=training_prs)
  mod1 <- glm(training_phenotype ~ ., data = df_train_1, family="binomial")
  
  # training AUC diff
  prediction <- prediction(predict(mod1, data=df_train_1, type="response"), training_phenotype)
  mod1_auc <- performance(prediction,measure="auc")@y.values[[1]]
  
  # for validation
  df_val_1 <- data.frame(prs=validation_prs)
  
  # validation AUC diff
  prediction <- prediction(predict(mod1, newdata=df_val_1, type="response"), validation_phenotype)
  mod1_val_auc <- performance(prediction,measure="auc")@y.values[[1]]
  
  x<-as.matrix(cbind(training_prs,training_theta_encoding))
  y<-training_phenotype
  n<-length(y)
  lam.seq <- exp(seq(log(lam.max),log(lam.min), length =nlambda))
  
  if(biglasso){ 
    x<-as.big.matrix(x)
    fit1 <- biglasso(x, y, family = "binomial",lambda=lam.seq)
    beta<-fit1$beta[-1,]
  }else{
    fit1 <- glmnet(x, y, family = "binomial",lambda=lam.seq)
    beta<-fit1$beta
  }
  
  training_validation_t<-matrix(0, nlambda, 4)
  support.list<-list()
  
  for(i in 1:nlambda){
    #if (i%%10==0){
    #  print(paste("nlambda",i))
    #}
    
    Theta<-beta[-1,i]
    support<-which(Theta!=0)
    
    if(length(support)>0){
      
      Theta_local<-numeric(p)
      y<-validation_phenotype
      x<-as.matrix(cbind(validation_prs,validation_theta_encoding[,support]))
      if(is.null(delta)){
        fit0<-glm(y~x, family="binomial")
        support<-support[!is.na(fit0$coefficients[-c(1:2)])]
        Theta_local[support]<-as.vector(summary(fit0)$coefficients[-c(1:2),3])
      }else{
        fit0 <- glmnet(x, y, family = "binomial",lambda =lam.seq[i])
        Theta1 <- as.vector(coef(fit0))
        fit0<-REF_DS_inf(x, y, family="binomial", lasso_est=Theta1, delta)
        Theta1<-as.vector(fit0$est/fit0$se)
        Theta_local[support]<-Theta1[-c(1:2)]
      }
      supp<-FDR_control(Theta, Theta_local,alpha)
      support1<-which(supp==1)
      support.list<-c(support.list, list(support1) )
      
      if(length(support1)>0){
        df_train_2 <- data.frame(prs=training_prs, theta_new = training_theta_encoding[,support1])
        df_val_2 <- data.frame(prs=validation_prs, theta_new = validation_theta_encoding[,support1])
        mod2 <- glm(training_phenotype ~ ., data = df_train_2, family="binomial")
        
        # training AUC diff
        prediction <- prediction(predict(mod2, data=df_train_2, type="response"), training_phenotype)
        mod2_auc <- performance(prediction,measure="auc")@y.values[[1]]
        training_validation_t[i,1] <- mod2_auc - mod1_auc
        
        # validation AUC diff
        prediction <- prediction(predict(mod2, newdata=df_val_2, type="response"), validation_phenotype)
        mod2_val_auc <- performance(prediction,measure="auc")@y.values[[1]]
        training_validation_t[i,2] <- mod2_val_auc - mod1_val_auc
        
        diff <- abs(training_validation_t[i,1] - training_validation_t[i,2])
        noise <- rnorm(1, 0, tolerance)
        abovethr <- (diff > threshold + noise)
        if (abovethr == FALSE){
          training_validation_t[i,2] <- training_validation_t[i,1]
        }else{
          noise <- rnorm(1, 0, tolerance)
          training_validation_t[i,2] <- training_validation_t[i,2] + noise
        }
        
        # selection of scores greater than err
        training_validation_t[i,3] <- (training_validation_t[i,1] > err) &
          (training_validation_t[i,2] > err)
        
        if( training_validation_t[i,3]==1){
          mod1_auc<-mod2_auc
          mod1_val_auc<-mod2_val_auc
        }
        training_validation_t[i,4]<- mod2_val_auc
      }
    }else{
      support.list<-c(support.list, list(0))
    }
  }
  
  index<-which(training_validation_t[,3]==1)
  if(length(index)==0){
    warning("error control is to large")
    index<-which.max(training_validation_t[,4])
  }
  
  a<-which.max(training_validation_t[index,4])
  support.F<-support.list[[index[a]]]
  
  return(support.F)
}

#' Adaptive Variable Selection for Continuous Outcomes
#'
#' This function identifies important \eqn{\theta_{SNPs}} through adaptive variable
#' selection and FDR control for continuous outcomes. 
#'
#' @param training_prs dataframe containing the training polygenic risk scores
#' @param validation_prs dataframe containing the validation polygenic risk scores
#' @param training_phenotype dataframe containing the training phenotype data
#' @param validation_phenotype dataframe containing the validation phenotype data
#' @param training_theta_encoding dataframe containing the training thetas
#' @param validation_theta_encoding dataframe containing the validation thetas
#' @param biglasso logical flag indicating whether to use big data optimization. 
#' If true, then the program uses \link[biglasso]{biglasso}. Or else, it uses 
#' the standard \link[glmnet]{glmnet} (default: FALSE)
#' @param lam.max maximum value of the regularization parameter (lambda) to be 
#' considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param lam.min minimum value of the regularization parameter (lambda) to be 
#' considered. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param nlambda number of different lambda values to be evaluated between 
#' lam.max and lam.min. This value is used to generate the lambda sequence in \link[glmnet]{glmnet}. 
#' @param alpha desired FDR control level.
#' @param tolerance tolerance level for noise in mirror statistics.
#' @param threshold threshold value for determining significance in mirror 
#' statistics.
#' @param err error threshold for model performance improvement.
#' @param delta additional tuning parameter for the FDR control.
#' @return A vector of selected \eqn{\theta_{SNPs}}. 
#' @export
adaptive_variable_selection_continuous<-function(training_prs, validation_prs, 
                                                 training_phenotype,validation_phenotype, 
                                                 training_theta_encoding, validation_theta_encoding, 
                                                 biglasso=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                                                 alpha=0.1, tolerance=0.025,threshold=0.01,err=1e-5, delta=NULL){
  
  p<-ncol(training_theta_encoding)
  
  # for training
  df_train_1 <- data.frame(prs=training_prs)
  mod1 <- glm(training_phenotype ~ ., data = df_train_1, family="gaussian")
  
  # training AUC diff
  prediction <-(predict(mod1, data=df_train_1, type="response")- training_phenotype)^2
  mod1_auc <- mean(prediction)
  
  
  # for validation
  df_val_1 <- data.frame(prs=validation_prs)
  
  # validation AUC diff
  prediction <- (predict(mod1, newdata=df_val_1, type="response")- validation_phenotype)^2
  mod1_val_auc <- mean(prediction)
  
  
  x<-as.matrix(cbind(training_prs,training_theta_encoding))
  y<-training_phenotype
  n<-length(y)
  #Generate a sequence of lambda values that are exponentially spaced between 
  #lam.max and lam.min with nlambda evenly spaced points in between. 
  lam.seq <- exp(seq(log(lam.max),log(lam.min), length =nlambda))
  
  if(biglasso){ 
    x<-as.big.matrix(x)
    fit1 <- biglasso(x, y, family = "gaussian",lambda=lam.seq)
    beta<-fit1$beta[-1,]
  }else{
    fit1 <- glmnet(x, y, family = "gaussian",lambda=lam.seq)
    beta<-fit1$beta
  }
  
  nlambda<- length(fit1$lambda)
  training_validation_t<-matrix(0, nlambda, 4)
  support.list<-list()
  
  for(i in 1:nlambda){
    #if (i%%10==0){
    #  print(paste("nlambda",i))
    #}
    
    Theta<-beta[-1,i]
    support<-which(Theta!=0)
    
    if(length(support)>0){
      
      Theta_local<-numeric(p)
      y<-validation_phenotype
      x<-as.matrix(cbind(validation_prs,validation_theta_encoding[,support]))
      
      if(is.null(delta)){
        fit0<-glm(y ~ x, family="gaussian")
        support<-support[!is.na(fit0$coefficients[-c(1:2)])]
        Theta_local[support]<-as.vector(summary(fit0)$coefficients[-c(1:2),3])
      }else{
        fit0 <- glmnet(x, y, family = "gaussian",lambda =lam.seq[i])
        Theta1 <- as.vector(coef(fit0))
        fit0<-REF_DS_inf(x, y, family="gaussian", lasso_est=Theta1, delta)
        Theta1<-as.vector(fit0$est/fit0$se)
        Theta_local[support]<-Theta1[-c(1:2)]
      }
      supp<-FDR_control(Theta, Theta_local,alpha)
      
      support1<-which(supp==1)
      support.list<-c(support.list, list(support1) )
      
      if(length(support1)>0){
        df_train_2 <- data.frame(prs=training_prs, theta_new = training_theta_encoding[,support1])
        df_val_2 <- data.frame(prs=validation_prs, theta_new = validation_theta_encoding[,support1])
        mod2 <- glm(training_phenotype ~ ., data = df_train_2, family="gaussian")
        
        # training AUC diff
        prediction <- (predict(mod2, data=df_train_2, type="response")- training_phenotype)^2
        mod2_auc <- mean(prediction)
        training_validation_t[i,1] <-  mod1_auc-mod2_auc
        
        # validation AUC diff
        prediction <- (predict(mod2, newdata=df_val_2, type="response")- validation_phenotype)^2
        mod2_val_auc <- mean(prediction)
        training_validation_t[i,2] <-  mod1_val_auc-mod2_val_auc
        
        diff <- abs(training_validation_t[i,1] - training_validation_t[i,2])
        noise <- rnorm(1, 0, tolerance)
        abovethr <- (diff > threshold + noise)
        if (abovethr == FALSE){
          training_validation_t[i,2] <- training_validation_t[i,1]
        }else{
          noise <- rnorm(1, 0, tolerance)
          training_validation_t[i,2] <- training_validation_t[i,2] + noise
        }
        
        
        # selection
        training_validation_t[i,3] <- (training_validation_t[i,1] > err) &
          (training_validation_t[i,2] > err)
        
        if( training_validation_t[i,3]==1){
          mod1_auc<-mod2_auc
          mod1_val_auc<-mod2_val_auc
        }
        training_validation_t[i,4]<- mod2_val_auc
      }
    }else{
      support.list<-c(support.list, list(0))
    }
  }
  
  index<-which(training_validation_t[,3]==1)
  
  if(length(index)==0){
    warning("error control is to large")
    index<-which.min(training_validation_t[,4])
  }
  
  a<-which.min(training_validation_t[index,4])
  support.F<-support.list[[index[a]]]
  
  return(support.F)
}