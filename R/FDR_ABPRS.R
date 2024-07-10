library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

###########################################
##
## Part III.  FDR_ABPRS for GLM and LM
##
###########################################


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

#' FDR control via mirror statistics from GLM (binary phenotype)
#'
#' This function conducts FDR control through mirror statistics. 
#'
#' @param prs_train dataframe containing the training polygenic risk scores
#' @param prs_val dataframe containing the validation polygenic risk scores
#' @param dat_train dataframe containing the training phenotype data
#' @param dat_val dataframe containing the validation phenotype data
#' @param theta_train dataframe containing the training thetas
#' @param theta_val dataframe containing the validation thetas
#' @param bigdata logical flag indicating whether to use big data optimization 
#' (default: FALSE)
#' @param lam.max maximum value of the regularization parameter (lambda) to be 
#' considered
#' @param lam.min minimum value of the regularization parameter (lambda) to be 
#' considered
#' @param nlambda number of different lambda values to be evaluated between 
#' lam.max and lam.min
#' @param alpha desired FDR control level
#' @param tolerance tolerance level for noise in mirror statistics
#' @param threshold threshold value for determining significance in mirror 
#' statistics
#' @param err error threshold for model performance improvement
#' @param delta additional tuning parameter for the FDR control
#' @return A large list with the following components:
#'    - support.F: selected false support
#'    - support.T: selected true support based on positive performance metrics for 
#'    both validation and training scores (removed later)
#'    - support.list: list of selected support for each lambda value.
#'    - training_validation_t: matrix of training and validation performance metrics (AUC score)
#' @export
FDR_selection_GLM<-function(prs_train=prs_train, prs_val=prs_val, dat_train=dat_train,
                            dat_val=dat_val, theta_train=theta_train, theta_val=theta_val,
                            bigdata=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                            alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-6, delta=NULL){
  
  p<-dim(theta_train)[2]
  
  # for training
  df_train_1 <- data.frame(prs=prs_train)
  mod1 <- glm(dat_train ~ ., data = df_train_1, family="binomial")
  
  # training AUC diff
  prediction <- prediction(predict(mod1, data=df_train_1, type="response"), dat_train)
  mod1_auc <- performance(prediction,measure="auc")@y.values[[1]]
  
  # for validation
  df_val_1 <- data.frame(prs=prs_val)
  
  # validation AUC diff
  prediction <- prediction(predict(mod1, newdata=df_val_1, type="response"), dat_val)
  mod1_val_auc <- performance(prediction,measure="auc")@y.values[[1]]
  
  x<-as.matrix(cbind(prs_train,theta_train))
  y<-dat_train
  n<-length(y)
  lam.seq <- exp(seq(log(lam.max),log(lam.min), length =nlambda))
  
  if(bigdata){ 
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
    if (i%%10==0){
      print(paste("nlambda",i))
    }
    
    Theta<-beta[-1,i]
    support<-which(Theta!=0)
    
    if(length(support)>0){
      
      Theta_local<-numeric(p)
      y<-dat_val
      x<-as.matrix(cbind(prs_val,theta_val[,support]))
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
        df_train_2 <- data.frame(prs=prs_train, theta_new = theta_train[,support1])
        df_val_2 <- data.frame(prs=prs_val, theta_new = theta_val[,support1])
        mod2 <- glm(dat_train ~ ., data = df_train_2, family="binomial")
        
        # training AUC diff
        prediction <- prediction(predict(mod2, data=df_train_2, type="response"), dat_train)
        mod2_auc <- performance(prediction,measure="auc")@y.values[[1]]
        training_validation_t[i,1] <- mod2_auc - mod1_auc
        
        # validation AUC diff
        prediction <- prediction(predict(mod2, newdata=df_val_2, type="response"), dat_val)
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
  
  a<-which.max(training_validation_t[,4])
  support.T<-support.list[[a]]
  
  result<-list(support.F=support.F,
               support.T=support.T,
               support.list=support.list,
               training_validation_t=training_validation_t )
  
  return(result)
}

#' FDR control via mirror statistics from LM (continuous phenotype)
#'
#' This function conducts FDR control through mirror statistics. It uses 
#' family = "gaussian" for the linear model. Data is optimized based on 
#' Mean Squared Error (MSE) reduction. 
#'
#' @param prs_train dataframe containing the training polygenic risk scores
#' @param prs_val dataframe containing the validation polygenic risk scores
#' @param dat_train dataframe containing the training phenotype data
#' @param dat_val dataframe containing the validation phenotype data
#' @param theta_train dataframe containing the training thetas
#' @param theta_val dataframe containing the validation thetas
#' @param bigdata logical flag indicating whether to use big data optimization. 
#' If true, then the program uses \link[biglasso]{biglasso}. Or else, it uses 
#' the standard \link[glmnet]{glmnet} (default: FALSE)
#' @param lam.max maximum value of the regularization parameter (lambda) to be 
#' considered
#' @param lam.min minimum value of the regularization parameter (lambda) to be 
#' considered
#' @param nlambda number of different lambda values to be evaluated between 
#' lam.max and lam.min
#' @param alpha desired FDR control level
#' @param tolerance tolerance level for noise in mirror statistics
#' @param threshold threshold value for determining significance in mirror 
#' statistics
#' @param err error threshold for model performance improvement
#' @param delta additional tuning parameter for the FDR control
#' @return a large list with the following components:
#'    - support.F: list of selected False SNPs
#'    - support.T: list of selected True SNPs
#'    - support.list: list of all selected SNPs
#'    - training_validation_t: list of AUC scores
#' @export
FDR_selection_LM<-function(prs_train=prs_train, prs_val=prs_val, dat_train=dat_train,
                           dat_val=dat_val, theta_train=theta_train, theta_val=theta_val, 
                           bigdata=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                           alpha=0.1, tolerance=0.025,threshold=0.01,err=1e-6, delta=NULL){
  
  p<-ncol(theta_train)
  
  # for training
  df_train_1 <- data.frame(prs=prs_train)
  mod1 <- glm(dat_train ~ ., data = df_train_1, family="gaussian")
  
  # training AUC diff
  prediction <-(predict(mod1, data=df_train_1, type="response")- dat_train)^2
  mod1_auc <- mean(prediction)
  
  
  # for validation
  df_val_1 <- data.frame(prs=prs_val)
  
  # validation AUC diff
  prediction <- (predict(mod1, newdata=df_val_1, type="response")- dat_val)^2
  mod1_val_auc <- mean(prediction)
  
  
  x<-as.matrix(cbind(prs_train,theta_train))
  y<-dat_train
  n<-length(y)
  #Generate a sequence of lambda values that are exponentially spaced between 
  #lam.max and lam.min with nlambda evenly spaced points in between. 
  lam.seq <- exp(seq(log(lam.max),log(lam.min), length =nlambda))
  
  if(bigdata){ 
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
    if (i%%10==0){
      print(paste("nlambda",i))
    }
    
    Theta<-beta[-1,i]
    support<-which(Theta!=0)
    
    if(length(support)>0){
      
      Theta_local<-numeric(p)
      y<-dat_val
      x<-as.matrix(cbind(prs_val,theta_val[,support]))
      
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
        df_train_2 <- data.frame(prs=prs_train, theta_new = theta_train[,support1])
        df_val_2 <- data.frame(prs=prs_val, theta_new = theta_val[,support1])
        mod2 <- glm(dat_train ~ ., data = df_train_2, family="gaussian")
        
        # training AUC diff
        prediction <- (predict(mod2, data=df_train_2, type="response")- dat_train)^2
        mod2_auc <- mean(prediction)
        training_validation_t[i,1] <-  mod1_auc-mod2_auc
        
        # validation AUC diff
        prediction <- (predict(mod2, newdata=df_val_2, type="response")- dat_val)^2
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
  
  a<-which.min(training_validation_t[,4])
  support.T<-support.list[[a]]
  
  result<-list(support.F=support.F,
               support.T=support.T, 
               support.list=support.list,
               training_validation_t=training_validation_t )
  
  return(result)
}