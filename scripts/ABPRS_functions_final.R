library(dplyr)
library(ROCR)
library(glmnet)
#library(biglasso)
#library(bigsnpr)
library(readr)
library(data.table)
library(Matrix)

###########################################
##
## Part I. Learning betas and covert theta
##
###########################################


load_geno <- function(dir_geno){
  geno <- fread(dir_geno,header=T)
  snp_names <- colnames(geno)
  snp_names <- snp_names[7:length(snp_names)]
  geno <- geno[,-c(1,3,4,5,6)]
  geno <- as.matrix(geno)
  ## zero imputation
  geno[is.na(geno)] <- 0
  
  colnames(geno) <- c("IID", snp_names)
  #print(dim(geno))
  snp_list <- as.data.frame(colnames(geno))
  snp_list <- as.data.frame(gsub("_.*","",snp_list[,1]))
  colnames(geno) <- snp_list[,1]
  geno <- as.matrix(geno)
  return(geno)
}


###########functions#########
learning_betas <- function(pheno=pheno_train, prs=prs_train, geno=snp_train, cov=cov_train, family="binomial"){
  p <- ncol(geno)
  betas <- data.frame(SNP=colnames(geno),beta1=rep(NA,p), beta2=rep(NA,p), pval=rep(NA,p))
  geno<-as.matrix(geno)
  for (i in 1:p){
    geno_f <- as.factor(geno[,i])
    if(length(levels(geno_f))>1){
      df <- data.frame(y = pheno, geno_f, prs)
      df <- cbind(df, cov)
      mod <- glm(y ~ ., data=df, family = family)
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


convert_theta <- function(geno, betas_mat){
  # for (i in 1:ncol(geno)){
  #   geno[,i] <- ifelse(geno[,i] == 0, 0,
  #                      ifelse(geno[,i] == 1, betas_mat[i,2], betas_mat[i,3]))
  # }
  geno<-t(t(I(geno==1))*betas_mat$beta1+t(I(geno==2))*betas_mat$beta2)
  return(geno)
}


###################################################
##
## Part II. AUC test for GLM and predict loss for LM 
##
###################################################

#######AUC on test data functions####
expit<-function(x){
  1/(1+exp(-x))
}

predict_test<-function(x,est){
  mu<-est[1]+x%*%est[-1]
  return(expit(mu))
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



AUC_test<-function(support,data_list, lasso=FALSE, lambda=NULL,delta=0){
  prs_train <-data_list$prs_train
  prs_test <- data_list$prs_test
  
  
  dat_train <- data_list$dat_train
  
  dat_test <-data_list$dat_test
  if(length(support)>0){
    df_train_2 <- data.frame(prs=prs_train, theta_new = data_list$theta_train[,support])
    df_test_2 <- data.frame(prs=prs_test, theta_new = data_list$theta_test[,support])
  }else{
    df_train_2 <- data.frame(prs=prs_train)
    df_test_2 <- data.frame(prs=prs_test)
  }
  
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
    mod2 <- glm(dat_train ~ ., data = df_train_2, family="binomial")
    # training AUC diff
    prediction <- prediction(predict(mod2, newdata=df_test_2, type="response"), dat_test)
    mod2_val_auc <- performance(prediction, measure="auc")@y.values[[1]]
  }
  return(mod2_val_auc)
  
}

####predict loss for LM on test data#####
Loss_test<-function(support,data_list,lasso=FALSE, lambda=NULL, delta=0){
  prs_train <-data_list$prs_train
  prs_test <-data_list$prs_test
  
  dat_train <- data_list$dat_train
  dat_test <-data_list$dat_test
  if(length(support)>0){
    df_train_2 <- data.frame(prs=prs_train, theta_new = data_list$theta_train[,support])
    df_test_2 <- data.frame(prs=prs_test, theta_new = data_list$theta_test[,support])
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


###########################################
##
## Part III.  FDR_ABPRS for GLM and LM
##
###########################################


FDR_control<-function(beta1, beta2, q){
  M<- sign(beta1*beta2)*(abs(beta1)+abs(beta2)) 
  tau.seq<-abs(M[I(M!=0)])
  
  tau.seq<-tau.seq[order(tau.seq)]
  for (tau in c(0,tau.seq)) {
    fdp<- sum(I(M<=(-tau)))/max(sum(I(M>=tau)),1)
    if(fdp<=q){
      tau.hat<-tau
      break
    } 
    tau.hat<-tau
  }
  
  p<-length(M)
  index<-which(M>=tau.hat)
  S<-rep(0,p)
  S[index]<-1
  return(S)
}

#####FDR-ABPRS functions for GLM. ###

FDR_selection<-function(data_list, bigdata=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                        alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-6, delta=NULL){
  
  prs_train <-data_list$prs_train
  prs_val <-data_list$prs_val
  
  theta_train <- data_list$theta_train
  theta_val <- data_list$theta_val
  
  dat_train <-data_list$dat_train
  dat_val <-data_list$dat_val
  
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
      support.list<-c(support.list, list(NULL))
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

#####FDR-ABPRS functions form LM. ###

FDR_selection_LM<-function(data_list, bigdata=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                           alpha=0.1, tolerance=0.025,threshold=0.01,err=1e-6, delta=NULL,sparse=FALSE){
  
  prs_train <- data_list$prs_train
  prs_val <- data_list$prs_val
  
  theta_train <- data_list$theta_train
  theta_val <- data_list$theta_val
  
  dat_train <- data_list$dat_train
  dat_val <- data_list$dat_val
  
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
  
  
  if(sparse){
    x<-Matrix(cbind(prs_train,theta_train),sparse = TRUE)
  }else{
    x<-as.matrix(cbind(prs_train,theta_train))
  }
  
  y<-dat_train
  n<-length(y)
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
      support.list<-c(support.list, list(NULL))
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




