library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)

#' Generate Simulation Data
#' 
#' This function generates simulated phenotype and genotype information.
#' 
#' @param m total number of snps
#' @param effect_snps_vec vector of length 4 containing the number of effective snps
#' with additive, dominant, recessive, and codominant encoding respectively
#' @param n total number of samples
#' @param effect_size effect size of the effective snps, which acts differently
#' depending on binary or continuous outcomes. See \eqn{\beta} in the equations sections. 
#' @param beta0 a value to ensure a certain percentage of case phenotype for 
#' binary outcomes and to serve as a baseline measurement for continuous outcomes. 
#' See \eqn{\beta_{0}} in the equations section. 
#' @param binary a boolean flag indicating whether the data has binary outcomes (TRUE)
#' or continuous outcomes (FALSE). 
#' @return A dataframe containing the simulated phenotype in the first column and the
#' simulated genotype encoding in the rest of the columns. 
#' @section Equations:
#' 
#' For binary outcomes, the phenotypes are simulated using a logistic model, as shown in the equation below: \cr
#' \deqn{\log \left(\frac{P(Y=1| \mathbf{X})}{1-P(Y=1| \mathbf{X})} \right)=\beta_0+\left(\mathbf{X}_{ADD}, \mathbf{X}_{Non-ADD} \right)\bm{\beta}}
#' 
#' For continuous outcomes, the phenotypes are simulated using the equation below: \cr
#' \deqn{{Y}=\beta_0+\left(\mathbf{X}_{ADD}, \mathbf{X}_{Non-ADD} \right)\bm{\beta}}
#' 
#' @export
data_simulation <- function(m, effect_snps_vec, n, maf, effect_size, beta0, binary){
  
  #Total number of snps
  effect_snps <- sum(effect_snps_vec)
  #Total number of non effective snps
  non_effect_snps <- m - effect_snps
  
  # Get effect SNPs
  dat_effect <- c()
  dat_effect_trans <- c()
  
  # Number of effective additive, dominant, recessive, codominant snps
  effect_add <- effect_snps_vec[1]
  effect_dom <- effect_snps_vec[2]
  effect_rec <- effect_snps_vec[3]
  effect_cod <- effect_snps_vec[4] 
  
  # generate ADD SNPs effect
  if  (effect_add > 0){
    for (i in 1:effect_add){
      x <- rbinom(n,1,maf)+rbinom(n,1,maf) #binomial distribution producing 0, 1, 2
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
      x[x==2] <- 0 #change 2 to 0 so that 0&2 (AA & aa) have no effect; Aa has 1 effect
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
  
  # Combined matrix with empty phenotype
  pheno <- matrix(nrow=n, ncol=1, 0)
  dat <- cbind(pheno, dat_effect, dat_non_effect)
  dat_trans <- cbind(pheno, dat_effect_trans, dat_non_effect)
  
  #Phenotype Generation
  if(binary){
    coef<-rnorm(effect_snps,effect_size,sd=0.01)
    linear_predictor = beta0 + dat_trans[, c(2:(effect_snps + 1))] %*% as.matrix(coef)
    prob <- exp(linear_predictor) / (1 + exp(linear_predictor))
    dat[, 1] <- apply(prob, 2, function(x) rbinom(length(x), 1, x))
  } else{
    coef<-rnorm(effect_snps,effect_size,sd=0.01)
    dat[,1]<-beta0+dat_trans[,c(2:(effect_snps+1))] %*% as.matrix(coef)+rnorm(n)
  }
  
  colnames(dat) <- c("Phenotype", paste0("SNP", c(1:m)))
  dat <- as.data.frame(dat)
  
  return(dat)
}