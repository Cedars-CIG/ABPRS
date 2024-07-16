library(dplyr)
library(ROCR)
library(glmnet)

#' Genotype Generation
#' 
#' @param snps number of total snps
#' @param effect_snps_vec vector containing the number of effective additive, 
#' dominant, recessive, and codominant snps
#' @param n number of individuals
#' @param maf minor allele frequency
#' @return list consisting of:
#' - dat_effect: dataframe containing effective snps in the columns and individuals in the rows
#' - dat_effect_trans: similar with dat_effect
#' - dat_non_effect: dataframe containing non-effective snps in the columns and individuals in the rows
geno_generation <- function(snps, effect_snps_vec, n, maf){
  
  #Total number of snps
  effect_snps <- sum(effect_snps_vec)
  
  # Get effect SNPs
  dat_effect <- c()
  dat_effect_trans <- c()
  
  #Total number of non effective snps
  non_effect_snps <- snps - effect_snps
  
  # Number of effective additive, dominant, recessive, codominant snps
  effect_add <- effect_snps_vec[1]
  effect_dom <- effect_snps_vec[2]
  effect_rec <- effect_snps_vec[3]
  effect_cod <- effect_snps_vec[4] 
  
  # generate ADD SNPs effect
  if  (effect_add > 0){
    for (i in 1:effect_add){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf) #binary distribution producing 0, 1, 2
      dat_effect <- cbind(dat_effect,x)
      dat_effect_trans <- cbind(dat_effect_trans,x) 
    }
  }
  # generate DOM SNPs effect
  if  (effect_dom > 0){
    for (i in 1:effect_dom){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf) #same process as above, but...
      dat_effect <- cbind(dat_effect,x)
      x[x==2] <- 1 #change all the 2s into 1s since Aa has same effect as AA
      dat_effect_trans <- cbind(dat_effect_trans,x)
    }
  }
  # generate REC SNPs effect
  if  (effect_rec > 0){
    for (i in 1:effect_rec){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf)
      dat_effect <- cbind(dat_effect,x)
      x[x==1] <- 0
      x[x==2] <- 1 #change 2 to 1 and 1 to 0, since AA and Aa has no effect while aa has effect
      dat_effect_trans <- cbind(dat_effect_trans,x)
    }
  }
  # generate COD SNPs effect
  if  (effect_cod > 0){
    for (i in 1:effect_cod){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf)
      dat_effect <- cbind(dat_effect,x)
      x[x==2] <- 0 #change 2 to 0 so that 0&2 are AA & aa, which have no effect; Aa has 1 effect
      dat_effect_trans <- cbind(dat_effect_trans,x)
    }
  }
  # get non_effect SNPs
  dat_non_effect <- c() 
  if  (non_effect_snps > 0){
    for (i in 1:non_effect_snps){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf)
      dat_non_effect <- cbind(dat_non_effect,x)
    }
  }
  
  return(list(dat_effect, dat_effect_trans, dat_non_effect))
}

pheno_generation <-function(dat_trans, dat, effect_snps, coef){
  linear_predictor = beta0 + dat_trans[, c(2:(effect_snps + 1))] %*% as.matrix(coef)
  prob <- exp(linear_predictor) / (1 + exp(linear_predictor))
  dat[, 1] <- apply(prob, 2, function(x) rbinom(length(x), 1, x))
  return(dat)
}

prs_generation <-function(dat_train, dat_test, dat_val){
  snps=(dim(dat_train)[[2]]-1)
  # GWAS betas, on train 
  gaws <- matrix(nrow=snps,ncol=2)
  for (i in 1:snps){
    gaws[i,1]<-summary(glm(dat_train[,1]~dat_train[,i+1],family="binomial"))$coef[2,1]
    gaws[i,2]<-i
  }
  
  # PRS calculation using ALL SNPs
  prs_train <- dat_train[,(gaws[,2]+1)]%*%as.matrix(gaws[,1])
  prs_test <- dat_test[,(gaws[,2]+1)]%*%as.matrix(gaws[,1])
  prs_val <- dat_val[,(gaws[,2]+1)]%*%as.matrix(gaws[,1])
  
  return(list(prs_train=prs_train, prs_test=prs_test, prs_val=prs_val))
}

#Split Data into Training, Validation, and Test sets
data_generation <- function(snps, effect_snps_vec, n, maf, effect_size, beta0){
  
  #Total number of snps
  effect_snps <- sum(effect_snps_vec)
  
  # Generate and retrieve genotype matrices
  geno_list <- geno_generation(snps, effect_snps_vec, n, maf)
  dat_effect = geno_list[[1]]
  dat_effect_trans = geno_list[[2]]
  dat_non_effect = geno_list[[3]]
  
  # Combined matrix with empty phenotype
  pheno <- matrix(nrow=n, ncol=1, 0)
  dat <- cbind(pheno, dat_effect, dat_non_effect)
  dat_trans <- cbind(pheno, dat_effect_trans, dat_non_effect)
  
  #Split training, validation, and test data 
  #training(33.33%), validation(3*16.67%), test(16.67%) data split
  sample_id<-sample(1:5,size=n,replace=T,prob=c(1/3,1/6,1/6,1/6,1/6))
  training<-which(sample_id==1)
  validation_1<-which(sample_id==2)
  #validation_2<-which(sample_id==3)
  #validation_3<-which(sample_id==4)
  testing_id<-which(sample_id==5)
  
  dat_train<-dat[training,]
  dat_val<-dat[validation_1,]
  #dat_val_2<-dat[validation_2,]
  #dat_val_3<-dat[validation_3,]
  dat_test<-dat[testing_id,]
  
  dat_train_trans<-dat_trans[training,]
  dat_val_trans<-dat_trans[validation_1,]
  #dat_val_2_trans<-dat_trans[validation_2,]
  #dat_val_3_trans<-dat_trans[validation_3,]
  dat_test_trans<-dat_trans[testing_id,]
  
  #Phenotype for each sample
  coef<-rnorm(effect_snps,effect_size,sd=0.01) #generate effect_snps(number of effective snps) random variables of normal distribution with mean of effect_size and standard deviation of 1
  dat_train = pheno_generation(dat_train_trans, dat_train, effect_snps, coef)
  dat_val = pheno_generation(dat_val_trans, dat_val, effect_snps, coef)
  dat_test = pheno_generation(dat_test_trans, dat_test, effect_snps, coef)
  
  dat_list <- list(dat_train, dat_val, dat_test, 
                   dat_train_trans, dat_val_trans, dat_test_trans )
  return(dat_list)
}

#Sort SNPs based on AUC diff in descending form
sort_snps <- function(pheno_train, prs_train, theta_train){
  snps = dim(theta_train)[[2]]
  order_t <- matrix(nrow=snps,ncol=2)
      # col 1: mod2_auc - mod1_auc
      # col 2: SNP id
  # mod1: pheno ~ prs
  # mod2: pheno ~ prs + theta_i
  mod1 <- glm(pheno_train ~ prs_train,family="binomial")
  prediction <- prediction(predict(mod1, type="response"), pheno_train)
  mod1_auc <- performance(prediction,measure="auc")@y.values[[1]]
  p_value<-numeric(snps)
  for (i in 1:snps){
    mod2 <- glm(pheno_train ~ prs_train + theta_train[,i], family="binomial")
    p_value[i]<-summary(mod2)$coefficients[3,4]
    prediction <- prediction(predict(mod2, type="response"), pheno_train)
    mod2_auc <- performance(prediction,measure="auc")@y.values[[1]]
    order_t[i,1] <- mod2_auc - mod1_auc
    order_t[i,2] <- i
  }
  
  order_t_sorted <- order_t[order(-order_t[,1]),]
  snps_order <- order_t_sorted[,2]
  
  return(snps_order)
}

#' Generate simulation data list necessary for data analysis 
#' 
#' @param snps number of total snps
#' @param effect_snp_vec vector of length 4 containing the number of effective snps
#' with additive, dominant, recessive, and heterozygous encoding respectively
#' @param n total number of samples
#' @param effect_size effect size of the effective snps
#' @param beta0 beta0 value
#' @return A list consisting of:
#' - snps_order: snps that are ordered
#' - theta_train, theta_val, theta_test
#' - pheno_train, pheno_val, pheno_test
#' - prs_train, prs_val, prs_test
#' - geno_train, geno_val, geno_test
#' @export
simulation_generation <- function(snps, effect_snps_vec, n, maf, effect_size, beta0){
  dat_list = data_generation(snps, effect_snps_vec, n, maf, effect_size, beta0)
  pheno_train = dat_list[[1]][,1]
  pheno_val = dat_list[[2]][,1]
  pheno_test = dat_list[[3]][,1]
  
  geno_train = dat_list[[1]][,-1]
  geno_val = dat_list[[2]][,-1]
  geno_test = dat_list[[3]][,-1]

  #PRS Generation
  prs_list = prs_generation(dat_list[[1]], dat_list[[2]], dat_list[[3]])
  prs_train = prs_list[[1]]
  prs_val = prs_list[[2]]
  prs_test = prs_list[[3]]
  
  #Generate covariates
  cov_train <- matrix(nrow=length(pheno_train), ncol=1, 0)
  
  #Learning betas: using prs_train and dat_trian
  betas_mat <- learning_betas(pheno_train, geno_train, prs_train, cov_train)
  
  #Convert thetas
  theta_train <-convert_theta(geno=geno_train, betas=betas_mat)
  theta_val <-convert_theta(geno=geno_val, betas=betas_mat)
  theta_test <- convert_theta(geno=geno_test, betas=betas_mat)
  
  #Order Snps based on AUC Diff
  snps_order=sort_snps(pheno_train, prs_train, theta_train)
  
  return(list(snps_order=snps_order, 
              theta_train=theta_train, theta_val=theta_val, theta_test=theta_test, 
              pheno_train=pheno_train, pheno_val=pheno_val, pheno_test=pheno_test, 
              prs_train=prs_train, prs_val=prs_val, prs_test=prs_test,
              geno_train=geno_train, geno_val=geno_val, geno_test=geno_test))
}
