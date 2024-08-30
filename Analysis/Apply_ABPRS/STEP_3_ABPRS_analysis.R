library(dplyr)
library(ROCR)

pheno <- "hypertension"

m_func <- function(data, race_num, df_coef){
  # filter race
  data <- data[data$race == race_num, ]
  
  y <- data$y
  X <- as.matrix(data[, -c(1,2,3)])
  
  df_X_names <- data.frame(var = colnames(X))
  df_coef <- inner_join(df_coef, df_X_names, by="var")
  
  X_all <- X[, df_coef$var]
  X_prs <- X_all[,1]
  X_theta <- X_all[,-c(1)]
  
  m1 <- X_prs * df_coef$coef[1]
  m2 <- X_theta %*% df_coef$coef[-1]
  
  return(list(m1, m2))
}

metric_func <- function(data, race_num, m1, m2){
  
  # discrete metric: AUC
  # continuous metric: R-squared
  data <- data[data$race == race_num, ]
  data$m1 <- m1
  data$m2 <- m2
  
  if (pheno %in% c("AD", "breast_cancer", "hypertension", "t2d")){
    mod1 <- glm(y ~ prs, data = data, family="binomial")
    mod2 <- glm(y ~ m1 + m2, data = data, family="binomial")
    
    prediction <- prediction(predict(mod1, newdata=data, type="response"), data$y)
    mod1_metric <- performance(prediction,measure="auc")@y.values[[1]]
    
    prediction <- prediction(predict(mod2, newdata=data, type="response"), data$y)
    mod2_metric <- performance(prediction,measure="auc")@y.values[[1]]
  } else {
    mod1 <- glm(y ~ prs, data = data, family="gaussian")
    mod2 <- glm(y ~ m1 + m2, data = data, family="gaussian")
    
    pred1 <- predict(mod1, newdata=data)
    mod1_metric <- cor(data$y, pred1)^2
    
    pred2 <- predict(mod2, newdata=data)
    mod2_metric <- cor(data$y, pred2)^2
  }
  df_rslt <- data.frame(external_PRS_metric=mod1_metric,
                        AB_PRS_metric=mod2_metric,
                        increase_percent=(mod2_metric-mod1_metric)/mod1_metric)
  return(df_rslt)
}

## loading data

prs_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/external_PRS/",pheno,"_external_PRS.csv")
data <- read.csv(prs_dir, header=T)


for (i in 1:22){
  #(paste("chr",i,"start..."))
  theta_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/theta_SNP/",pheno,"_chr",i,"_theta")
  
  if (file.exists(theta_dir)){
    theta <- read.table(theta_dir, header=T)
    data <- cbind(data, theta)
  }
}

pheno_dir <- "/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/AoU_phenotype/phenotype.csv"

if (pheno == "breast_cancer"){
  cols_select <- c("IID", "gender", "race", pheno)
  df_pheno <- read.csv(pheno_dir, header=T) %>% select(any_of(cols_select))
  colnames(df_pheno) <- c("IID", "gender", "race", "y")
  #print(dim(df_pheno))
  ## only select female
  df_pheno <- df_pheno[df_pheno$gender == 1,]
  #print(dim(df_pheno))
  df_pheno <- df_pheno %>% select(any_of(c("IID", "race", "y")))
} else{
  cols_select <- c("IID", "race", pheno)
  df_pheno <- read.csv(pheno_dir, header=T) %>% select(any_of(cols_select))
  colnames(df_pheno) <- c("IID", "race", "y")
}
#print("loading phenotype done.")

## filter continuous phenotype
if (pheno %in% c("bmi", "cholesterol", "hdl", "ldl")){
  dir <- "/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/AoU_phenotype/continuous_phenotype_range.csv"
  df_range <- read.csv(dir, header=T)
  df_range <- df_range[df_range$phenotype == pheno,]
  #print(df_range)
  # filter
  #print(paste("number of participants, before filter:", nrow(df_pheno)))
  df_pheno <- df_pheno[(df_pheno$y >= df_range$min)&(df_pheno$y <= df_range$max),]
  #print(paste("number of participants, after filter:", nrow(df_pheno)))
} else {
  #print("discrete phenotype!")
}

## inner_join phenotype and genotype
data <- inner_join(df_pheno, data, by="IID")

coef_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/UKBB_glm_model/",pheno,"_glm_model.RData")
load(coef_dir)

# if "Female":
#     gender = 1
# elif "Male":
#     gender = 2
# else:
#     gender = 3

# "Black or African American":
#     race = 1
# "White":
#     race = 2
# "Asian":
#     race = 3
# "Native Hawaiian or Other Pacific Islander":
#     race = 4
# "Middle Eastern or North African":
#     race = 5
#  others race 6

race_vec <- c("Black", "White", "Asian", "Hawaiian_Islander", "Mid_Eest_North_Africa", "other")
df_race <- data.frame(Race=race_vec)

rslt <- m_func(data, 1, coef_mod2)
df_rslt <- metric_func(data, 1, rslt[[1]], rslt[[2]])
for (i in 2:6){
  rslt <- m_func(data, i, coef_mod2)
  temp <- metric_func(data, i, rslt[[1]], rslt[[2]])
  df_rslt <- rbind(df_rslt, temp)
}
df_rslt <- cbind(df_race, df_rslt)

rslt_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/AB_PRS_result/",pheno,"_V3_result.csv")
write.csv(df_rslt, rslt_dir, row.names = FALSE)



