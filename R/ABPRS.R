#' ABPRS Function
#' 
#' This function runs the ABPRS method and returns the weights of 
#' the pre-trained prs and selected \eqn{\theta_{SNPs}} used for ABPRS calculation.
#' For the selected \eqn{\theta_{SNPs}}, it also returns their corresponding 
#' theta1 (theta encoding for heterozygous alleles Aa) and theta2 (theta encoding 
#' for homozygous alternative alleles aa) values, which are important for encoding
#' the genotype matrices when producing adaptively-boosted polygenic risk scores (AB-PRSs). 
#' 
#' @param pre_trained_prs dataframe containing the pre-trained training polygenic risk scores
#' @param validation_prs dataframe containing the pre-trained validation polygenic risk scores
#' @param training_phenotype vector containing the training phenotype data
#' @param validation_phenotype vector containing the validation phenotype data
#' @param training_genotype dataframe containing the training genotype data
#' @param validation_genotype dataframe containing the validation genotype data
#' @param family string specifying the type of model to fit a \link[stats]{glm}.
#' Use "binomial" for binary outcomes and "gaussian" for continuous outcomes. 
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
#' @param threshold threshold value ``v” for determining the significance of 
#' difference between training improvement and validation improvement, 
#' see equation (3) of ABPRS paper.
#' @param err error threshold for model performance improvement.
#' @param delta if ``NULL” then using GLM to fit the regression model on the validation dataset (default). 
#' If a small value is specified, the debiased lasso method is applied, where the Hessian matrix is adjusted 
#' by adding an identity matrix weighted by delta to ensure invertibility.
#' @return A dataframe with three columns. The first column named "Weight" stores 
#' the weights for pre-trained PRS scores and selected \eqn{\theta_{SNPs}}. The 
#' second and third columns stores the theta1 (theta encoding for heterozygous 
#' alleles Aa) and theta2 (theta encoding for homozygous alternative alleles aa)
#' of the selected \eqn{\theta_{SNPs}}, which are useful for encoding the genotype
#' when producing adaptively-boosted polygenic risk scores (AB-PRSs). 
#' @export
ABPRS <- function(pre_trained_prs, validation_prs, 
                   training_phenotype, validation_phenotype, 
                   training_genotype, validation_genotype,
                   family, biglasso=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                   alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-5, delta=NULL){
  
  if(family != "binomial" && family != "continuous"){
    warning("Different family than expected. Use \"binomial\" for binary outcomes and \"gaussian\" for continuous outcomes.")
  }
  
  # Learning beta coefficients using training data
  thetas_mat <- learning_theta_snps(phenotype=training_phenotype, genotype=training_genotype, 
                                    pre_trained_prs=pre_trained_prs, cov=NULL, family=family)
  
  # Convert thetas
  theta_train <-encoding_theta_snps(genotype=training_genotype, thetas=thetas_mat)
  theta_val <-encoding_theta_snps(genotype=validation_genotype, thetas=thetas_mat)

  support_ABPRS = adaptive_variable_selection(pre_trained_prs = pre_trained_prs, 
                                              validation_prs = validation_prs,
                                              training_phenotype = training_phenotype, 
                                              validation_phenotype = validation_phenotype,
                                              training_theta_encoding = theta_train, 
                                              validation_theta_encoding = theta_val,
                                              family=family, biglasso = biglasso, 
                                              lam.max = lam.max, lam.min = lam.min, nlambda = nlambda,
                                              alpha = alpha, tolerance = tolerance, threshold = threshold,
                                              err = err, delta = delta)
  
  # Generate model and Retrieve Weights
  trainABPRS <- data.frame(Pre_Trained_PRS=pre_trained_prs, theta_train[,support_ABPRS])
  modABPRS <- glm(training_phenotype ~ ., data = trainABPRS, family=family)
  
  #Return Dataframe with Weights and Betas
  select_thetas <- thetas_mat[thetas_mat$SNP %in% support_ABPRS,]
  weights <- data.frame(Weight=modABPRS$coefficients[-1], 
                        theta1=c(NA, select_thetas$theta1),
                        theta2=c(NA, select_thetas$theta2))
  return(weights)
}

#' Applying Weights to Produce AB-PRSs
#' 
#' This function applies the weights generated from the ABPRS function to produce
#' adaptively-boosted polygenic risk scores. 
#' 
#' @param pre_trained_prs a dataframe containing the pre-trained polygenic risk scores
#' @param genotype a dataframe containing the corresponding genotype information, 
#' whose rows represent individuals and columns represent SNPs. 
#' @param weights a dataframe containing the weights of the pre-trained prs and 
#' important \eqn{\theta_{SNPs}}, and the theta1 and theta2 values for those SNPs. 
#' This dataframe can be generated from \link[ABPRS]{ABPRS}. 
#' @return A dataframe with the adaptively-boosted polygenic risk scores. 
#' @export
apply_weights<- function(pre_trained_prs, genotype, weights){
  # Apply the model to generate polygenic risk scores for the data
  support_ABPRS <-rownames(weights)[-1]
  theta <- encoding_theta_snps(genotype[,support_ABPRS], weights[-1,])
  theta_mat <- data.frame(prs=pre_trained_prs, theta[, support_ABPRS])
  abprs = as.vector(as.matrix(theta_mat)%*%as.matrix(weights$Weight))
  abprs <- data.frame(ABPRS=abprs, row.names=rownames(genotype))
  return(abprs)
}
