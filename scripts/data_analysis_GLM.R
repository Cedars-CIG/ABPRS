source("/home/hujie/ABPRS/real-code/ABPRS_functions_all.R")

######hypertension####
# geno_name<-"hypertension"
# prs_name<-"hypertension_PGS002603"

####### t2d #######
# geno_name<-"t2d"
# prs_name<-"t2d_PGS000805"

####### AD #######
geno_name<-"AD_DNAnexus"
prs_name<-"AD_DNAnexus_PGS001775"

####### breast_cencer #######
# geno_name<-"breast_cancer"
# prs_name<-"breast_cancer_PGS000004"


# ####### chol #####
# geno_name<-"cholesterol"
# prs_name<-"cholesterol_PGS002286"
# 
# ####### ldl #####
# geno_name<-"ldl"
# prs_name<-"ldl_PGS000846"


#############  
prs_all<-fread(paste0("/home/hujie/ABPRS/prs/",prs_name,"_prs"))
pheno<-fread(paste0("/home/hujie/ABPRS/Phenotype/",geno_name,"_phenotype"), header = TRUE)

IID<-fread(paste0("/home/hujie/ABPRS/IID_complete/",geno_name,"_IID"))
IID_raw<-fread(paste0("/home/hujie/ABPRS/IID_complete/",geno_name,"_IID_raw"))
match_id <- match(IID_raw$V1,IID$V1)
rm(IID,IID_raw)
pheno<-pheno[match_id,]

train_id<-fread(paste0("/home/hujie/ABPRS/IID_split/",geno_name,"_train_id"))
val_id<-fread(paste0("/home/hujie/ABPRS/IID_split/",geno_name,"_val_id"))
test_id<-fread(paste0("/home/hujie/ABPRS/IID_split/",geno_name,"_test_id"))


##############


prs_train<-subset(prs_all, prs_all$V1%in% train_id$V1)
prs_val<-subset(prs_all, prs_all$V1%in% val_id$V1)
prs_test<-subset(prs_all, prs_all$V1%in% test_id$V1)
rm(prs_all)

dat_train<-subset(pheno,pheno$IID %in% train_id$V1)
dat_val<-subset(pheno,pheno$IID %in% val_id$V1)
dat_test<-subset(pheno,pheno$IID %in% test_id$V1)
rm(pheno)

dat_test<-dat_test$phenotype
dat_train<-dat_train$phenotype
dat_val<-dat_val$phenotype


prs_test<-prs_test$V2
prs_train<-prs_train$V2
prs_val<-prs_val$V2


#### obtain data ####

setwd(paste0("/home/hujie/ABPRS/data/AD"))


theta_T<-NULL
theta_F<-NULL

pval<-0.001
cat("pvalue",pval,"\n")
p1<-50
cat("lambda_min_dimension",p1,"\n")

### k stand for chr k #####
for (k in c(1:22)) {

name<-paste0(geno_name,"_chr",k,".txt")  
beta<-fread(name) 

file_name<-paste0(geno_name,"_chr",k,".raw")

ind<-which(beta$pval<pval)
# if(length(ind)>1500 ){
# ind<-order(beta$pval)[1:1500]
# }else 
# if(length(ind)<200){
#   ind<-order(beta$pval)[1:200]
# }
sel_set<-c(2,6+ind)
SNPs<-fread(file_name,select = sel_set, header = TRUE)

SNPs[is.na(SNPs)]<-0
IID<-SNPs$IID
SNPs<-SNPs[,-c("IID")]

betas_mat<-beta[ind,]
theta_snp<-convert_theta(as.matrix(SNPs),betas_mat)
rm(SNPs)

######## clean data


theta_train <- theta_snp[IID %in% train_id$V1,]
theta_test <- theta_snp[IID %in% test_id$V1,]
theta_val <- theta_snp[IID %in% val_id$V1,]



#########obtain data_list######
data_list<-NULL
data_list$prs_train<-prs_train
data_list$prs_val<-prs_val
data_list$prs_test<-prs_test

data_list$dat_train<-dat_train
data_list$dat_val<-dat_val
data_list$dat_test<-dat_test


data_list$theta_train<-theta_train
data_list$theta_val<-theta_val
data_list$theta_test<-theta_test


#############PRS only###
test_auc<-AUC_test(NULL,data_list)
test_auc
###lasso selection###
lam.max=3e-3
lam.min=2e-3
a<-5
while (a>4) {
  support_lasso<-Lasso_selection(data_list,family="binomial",lambda=lam.max)
  a<-length(support_lasso)
  print(a)
  if(a>4){
  lam.max<-2*lam.max
  }
}


if(ncol(theta_train)>2*p1){
  a<-0
  while (a<=p1) {
    support_lasso<-Lasso_selection(data_list,family="binomial",lambda=lam.min)
    a<-length(support_lasso)
    print(a)
    if(a<=p1){
      lam.min=lam.min/2
    }
  }
}else{
  p0<-ncol(theta_train)/2
  a<-0
  while (a<=p0) {
    support_lasso<-Lasso_selection(data_list,family="binomial",lambda=lam.min)
    a<-length(support_lasso)
    print(a)
    if(a<=p0){
      lam.min=lam.min/2
    }
  }
}
#Lasso_auc<-AUC_test(support_lasso,data_list,lasso=TRUE, lambda=lam.min,delta=1e-6)
#Lasso_auc
####FDR+ABPRS selection
t1<-Sys.time()
fit<-FDR_selection(data_list, bigdata=FALSE, lam.max=lam.max, lam.min=lam.min,nlambda=40,
                   alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-5, delta=NULL)
Sys.time()-t1

support_FDR.T<-fit$support.T
b<-colnames(as.matrix(theta_snp[,support_FDR.T,drop=FALSE]))

cat("chr",k, "length", length(support_FDR.T), "support", b, "\n")

FDR_auc<-AUC_test(support_FDR.T, data_list)
FDR_auc
cat("aucT",(FDR_auc-test_auc)/test_auc,"\n")

support_FDR.F<-fit$support.F

b<-colnames(theta_snp[,support_FDR.F,drop=FALSE])

cat("chr",k, "length", length(support_FDR.F), "support", b, "\n")

FDR_auc<-AUC_test(support_FDR.F, data_list)
FDR_auc
cat("aucF",(FDR_auc-test_auc)/test_auc,"\n")
######## save result of selected theta_snps

theta_T<-cbind(theta_T,theta_snp[,support_FDR.T,drop=FALSE])

theta_F<-cbind(theta_F,theta_snp[,support_FDR.F,drop=FALSE])

rm(data_list, theta_snp, theta_train,theta_test,theta_val)
}



###save the final result
name<-paste0(geno_name,"_theta_T_pval",pval,"p1",p1,".RData")
save(theta_T, file = name)



theta_train <- theta_T[IID %in% train_id$V1,]
theta_test <- theta_T[IID %in% test_id$V1,]
theta_val <- theta_T[IID %in% val_id$V1,]
data_list<-NULL
data_list$prs_train<-prs_train
data_list$prs_val<-prs_val
data_list$prs_test<-prs_test

data_list$dat_train<-dat_train
data_list$dat_val<-dat_val
data_list$dat_test<-dat_test


data_list$theta_train<-theta_train
data_list$theta_val<-theta_val
data_list$theta_test<-theta_test

support_FDR<-1:ncol(theta_train)
FDR_auc<-AUC_test(support_FDR,data_list)
cat("auc_all_T",(FDR_auc-test_auc)/test_auc,"\n")


name<-paste0(geno_name,"_theta_F_pval",pval,"p1",p1,".RData")
save(theta_F, file = name)

theta_train <- theta_F[IID %in% train_id$V1,]
theta_test <- theta_F[IID %in% test_id$V1,]
theta_val <- theta_F[IID %in% val_id$V1,]
data_list<-NULL
data_list$prs_train<-prs_train
data_list$prs_val<-prs_val
data_list$prs_test<-prs_test

data_list$dat_train<-dat_train
data_list$dat_val<-dat_val
data_list$dat_test<-dat_test


data_list$theta_train<-theta_train
data_list$theta_val<-theta_val
data_list$theta_test<-theta_test

support_FDR<-1:ncol(theta_train)
FDR_auc<-AUC_test(support_FDR,data_list)
cat("auc_all_F",(FDR_auc-test_auc)/test_auc,"\n")


save.image("workspace[1].RData")

############BHq p_value###

# t1<-Sys.time()
# support_BHq<-BH_selection(data_list,family="binomial", alpha = 0.1, lambda=3e-3)$selected_snps
# support_BHq
# Sys.time()-t1
# 
# AUC_test(support_BHq, data_list, lasso=TRUE, lambda=3e-3, delta=1e-4)
# 


