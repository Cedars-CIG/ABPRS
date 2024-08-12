#' ABPRS Function
#' 
#' This function adaptively-boosts pre-trained polygenic risk scores of a target
#' sample and returns the newly adaptively boosted polygenic risk scores of that
#' sample. 
#' 
#' @param target_prs dataframe containing the pre-trained prs that we wish to 
#' convert
#' @param target_genotype dataframe containing the genotype information of the 
#' target sample
#' @param training_prs dataframe containing the training polygenic risk scores
#' @param validation_prs dataframe containing the validation polygenic risk scores
#' @param training_phenotype dataframe containing the training phenotype data
#' @param validation_phenotype dataframe containing the validation phenotype data
#' @param training_genotype dataframe containing the training genotype data
#' @param validation_genotype dataframe containing the validation genotype data
#' @param family string specifying whether the data has binary or continuous outcomes.
#' Use "binary" for binary outcomes and "continuous" for continuous outcomes. 
#' (default: "binary"). 
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
#' @return Adaptively-boosted polygenic risk scores of the target sample.
#' @export
ABPRS <- function(target_prs, target_genotype,
                  training_genotype, validation_genotype,
                  training_prs, validation_prs, 
                  training_phenotype, validation_phenotype, 
                  family = "binary", 
                  biglasso=FALSE, lam.max=2e-3, lam.min=6e-5,nlambda=50,
                  alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-5, delta=NULL){
  
  if(family=="continuous"){
    family <- "gaussian"
  }else{
    family <- "binomial"
  }
  
  # Learning beta coefficients using training data
  thetas_mat <- learning_theta_snps(phenotype=training_phenotype, genotype=training_genotype, 
                                    prs=training_prs, cov=NULL, family=family)
  
  # Convert thetas
  theta_train <-encoding_theta_snps(genotype=training_genotype, thetas=thetas_mat)
  theta_val <-encoding_theta_snps(genotype=validation_genotype, thetas=thetas_mat)
  target_theta <- encoding_theta_snps(genotype=target_genotype, thetas=thetas_mat)
  
  if(family=="gaussian"){
    support_ABPRS = adaptive_variable_selection_continuous(training_prs = training_prs, 
                                                       validation_prs = validation_prs,
                                                       training_phenotype = training_phenotype, 
                                                       validation_phenotype = validation_phenotype,
                                                       training_theta_encoding = theta_train, 
                                                       validation_theta_encoding = theta_val,
                                                       biglasso = biglasso,
                                                       lam.max = lam.max, lam.min = lam.min, nlambda = nlambda,
                                                       alpha = alpha, tolerance = tolerance, threshold = threshold,
                                                       err = err, delta = delta)
  }else{
    support_ABPRS = adaptive_variable_selection_binary(training_prs = training_prs, 
                                                       validation_prs = validation_prs,
                                                       training_phenotype = training_phenotype, 
                                                       validation_phenotype = validation_phenotype,
                                                       training_theta_encoding = theta_train, 
                                                       validation_theta_encoding = theta_val,
                                                       biglasso = biglasso,
                                                       lam.max = lam.max, lam.min = lam.min, nlambda = nlambda,
                                                       alpha = alpha, tolerance = tolerance, threshold = threshold,
                                                       err = err, delta = delta)
  }
  
  # Generate binomial model from training data
  trainABPRS <- data.frame(prs=training_prs, theta_train[,support_ABPRS])
  modABPRS <- glm(training_phenotype ~ ., data = trainABPRS, family=family)
  
  # Apply the model to generate polygenic risk scores for testing data
  resultABPRS <- data.frame(prs=target_prs, target_theta[,support_ABPRS])
  ABPRS = as.vector(as.matrix(resultABPRS)%*%modABPRS$coefficients[-1])
  
  return(ABPRS)
}
