## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE, include = FALSE---------------------------------------
library(ABPRS)
library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)
library(ggplot2)
library(plotly)
library(tableHTML)
library(htmlwidgets)
library(DT)

## ----eval=FALSE---------------------------------------------------------------
# install.packages(“remotes”)
# remotes::install_github("Cedars-CIG/ABPRS", dependencies = TRUE)

## -----------------------------------------------------------------------------
library(ABPRS)

## ----eval=FALSE---------------------------------------------------------------
# # Core dependencies
# dplyr,
# ROCR,
# glmnet,
# biglasso,
# # Evaluation and visualization
# ggplot2,
# plotly,
# tableHTML,
# htmlwidgets

## ----data simulation----------------------------------------------------------
data <- data_simulation(m=250, effect_snps_vec = c(25, 0, 0, 25), n=10000, 
                        maf=0.3, effect_size=0.7, beta0= (-19.5), binary = TRUE)

## ----echo=FALSE---------------------------------------------------------------
str(data, list.len=0)
knitr::kable(data[1:8,1:11], row.names=TRUE, format="markdown")

## -----------------------------------------------------------------------------
mean(data[,1])

## ----sampling-----------------------------------------------------------------
sample <- sample(1:3,size=10000,replace=T,prob=c(0.5, 0.25, 0.25))
data_train = data[which(sample==1),]
data_val = data[which(sample==2),]
data_test = data[which(sample==3),]
rm(sample)

## -----------------------------------------------------------------------------
training_phenotype = data_train[,1]
validation_phenotype = data_val[,1]
testing_phenotype = data_test[,1]

training_genotype = data_train[,-1]
validation_genotype = data_val[,-1]
testing_genotype = data_test[,-1]

rm(data_test, data_train, data_val)

## -----------------------------------------------------------------------------
#Generate GWAS on training dataset
gwas <- matrix(nrow=250,ncol=2)
for (i in 1:250){
  gwas[i,1]<-summary(glm(training_phenotype~training_genotype[,i],family="binomial"))$coef[2,1]
  gwas[i,2]<-i
}
  
# PRS calculation
training_prs <- as.matrix(training_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])
validation_prs <- as.matrix(validation_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])
testing_prs <- as.matrix(testing_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])

## -----------------------------------------------------------------------------
str(training_prs, give.attr = FALSE)
str(validation_prs, give.attr = FALSE)
str(testing_prs, give.attr = FALSE)

## -----------------------------------------------------------------------------
weights <- ABPRS(pre_trained_prs = training_prs, validation_prs=validation_prs, 
                 training_genotype=training_genotype, validation_genotype=validation_genotype,
                 training_phenotype=training_phenotype, validation_phenotype=validation_phenotype,
                 family="binomial", covariate=NULL, biglasso=FALSE, 
                 lam.max=NULL, lam.min=NULL, nlambda=100,
                 alpha=0.1, tolerance=0.025, threshold=0.01,
                 err=1e-5, delta=NULL)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(head(weights[1:5, 1:3]), format="markdown")

## -----------------------------------------------------------------------------
rownames(weights)[-1]

## -----------------------------------------------------------------------------
testing_abprs <- apply_weights(testing_prs, testing_genotype, weights)

## -----------------------------------------------------------------------------
str(testing_abprs)

## -----------------------------------------------------------------------------
model_evaluation(phenotype=testing_phenotype, 
                 'Pre-Trained PRS' = testing_prs, 
                 'AB-PRS' = testing_abprs,
                 binary = TRUE, 
                 bin=10,
                 filename="Evaluation_Binary")

## -----------------------------------------------------------------------------
data <- data_simulation(m=250, effect_snps_vec = c(25, 0, 0, 25), n=10000, 
                        maf=0.3, effect_size=0.3, beta0= (-1), binary = FALSE)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(data[1:8,1:9], format="markdown")

## -----------------------------------------------------------------------------
sample <- sample(1:3,size=10000,replace=T,prob=c(0.5, 0.25, 0.25))
data_train = data[which(sample==1),]
data_val = data[which(sample==2),]
data_test = data[which(sample==3),]
rm(sample)

training_phenotype = data_train[,1]
validation_phenotype = data_val[,1]
testing_phenotype = data_test[,1]

training_genotype = data_train[,-1]
validation_genotype = data_val[,-1]
testing_genotype = data_test[,-1]

rm(data_test, data_train, data_val)

## -----------------------------------------------------------------------------
#Generate GWAS on training dataset
gwas <- matrix(nrow=250,ncol=2)
for (i in 1:250){
  gwas[i,1]<-summary(glm(training_phenotype~training_genotype[,i],family="gaussian"))$coef[2,1]
  gwas[i,2]<-i
}
  
# PRS calculation
training_prs <- as.matrix(training_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])
validation_prs <- as.matrix(validation_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])
testing_prs <- as.matrix(testing_genotype[,(gwas[,2])])%*%as.matrix(gwas[,1])

## -----------------------------------------------------------------------------
thetas_mat <- learning_theta_snps(phenotype=training_phenotype, genotype=training_genotype,
                                  pre_trained_prs=training_prs, covariate=NULL, family="gaussian")

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(thetas_mat[1:8,1:4], format="markdown")

## -----------------------------------------------------------------------------
theta_train <- encoding_theta_snps(training_genotype, thetas_mat)
theta_val <- encoding_theta_snps(validation_genotype, thetas_mat)
theta_test <- encoding_theta_snps(testing_genotype, thetas_mat)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(theta_train[1:8,1:8], format="markdown")

## -----------------------------------------------------------------------------
support_ABPRS <- adaptive_variable_selection(pre_trained_prs = training_prs, 
                                            validation_prs = validation_prs,
                                            training_phenotype = training_phenotype, 
                                            validation_phenotype = validation_phenotype,
                                            training_theta_encoding = theta_train, 
                                            validation_theta_encoding = theta_val,
                                            family = "gaussian", biglasso = FALSE,
                                            lam.max = NULL, lam.min = NULL, nlambda = 100,
                                            alpha = 0.1, tolerance = 0.025, threshold = 0.01,
                                            err = 1e-05, delta = NULL)

## -----------------------------------------------------------------------------
support_ABPRS

## -----------------------------------------------------------------------------
# Generate gaussian model from training data
trainABPRS <- data.frame(prs=training_prs, theta_train[,support_ABPRS])
modABPRS <- glm(training_phenotype ~ ., data = trainABPRS, family="gaussian")

## -----------------------------------------------------------------------------
# Apply the model to generate polygenic risk scores for testing data
testABPRS <- data.frame(prs=testing_prs, theta_test[,support_ABPRS])
testing_abprs = as.matrix(testABPRS)%*%modABPRS$coefficients[-1]

## -----------------------------------------------------------------------------
str(as.vector(testing_abprs))

## -----------------------------------------------------------------------------
model_evaluation(phenotype=testing_phenotype, 
                 'Pre-Trained PRS' = testing_prs, 
                 'AB-PRS' = testing_abprs,
                 binary = FALSE, 
                 bin=10,
                 filename="Evaluation_Continuous")

