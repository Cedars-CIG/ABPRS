source("/home/hujie/ABPRS/real-code/FDR_AB_PRS_functions[2].R")

args <- commandArgs(trailingOnly = TRUE)
print(args)
####### chol #####
#geno_name<-"cholesterol"
#prs_name<-"cholesterol_PGS002286"

# ###### ldl #####
# geno_name<-"ldl"
# prs_name<-"ldl_PGS000846"

# # ###### hdl #####
# geno_name<-"hdl"
# prs_name<-"hdl_PGS000845"

# # ###### bmi #####
 geno_name<-"bmi"
 prs_name<-"bmi_PGS000841"

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
setwd("/home/hujie/ABPRS/data/bmi")


theta_T<-NULL
theta_F<-NULL
pval<-as.numeric(args[1])
cat("p_value",pval,"\n")
p1<-as.numeric(args[2])
cat("lambda_min_dimension",p1,"\n")

p2<-as.numeric(args[3])
cat("max_dimension",p2,"\n")

### k stand for chr k #####
for (k in c(1:22)) {

name<-paste0(geno_name,"_chr",k,".txt")  
beta<-fread(name) 
file_name<-paste0(geno_name,"_chr",k,".raw")
ind<-which(beta$pval<pval)
if(length(ind)>p2 ){
ind<-order(beta$pval)[1:p2]
}

dim_snp<-fread(file_name,nrows = 1,header = TRUE)
ind1<-which(colnames(dim_snp)[-c(1:6)] %in%  beta$SNP[ind])

#else 
#if(length(ind)<100){
#  ind<-order(beta$pval)[1:100]
#}
sel_set<-c(2,6+ind1)
SNPs<-fread(file_name,select = sel_set, header = TRUE)
SNPs[is.na(SNPs)]<-0
IID<-SNPs$IID
SNPs<-SNPs[,-c("IID")]

betas_mat<-beta[ind,]
theta_snp<-convert_theta(as.matrix(SNPs),betas_mat)
rm(SNPs)
gc()
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

rm(theta_train,theta_test,theta_val)
gc() 
#############PRS only###
test_loss<-Loss_test(NULL,data_list)
test_loss
###lasso selection###
lam.max=6e-3
lam.min=2e-3
a<-5
while (a>4) {
  support_lasso<-Lasso_selection(data_list,family="gaussian",lambda=lam.max)
  a<-length(support_lasso)
  print(a)
  if(a>4){
  lam.max<-2*lam.max
  }
}


if(ncol(theta_snp)>p1){
a<-0
while (a<=p1) {
support_lasso<-Lasso_selection(data_list,family="gaussian",lambda=lam.min)
a<-length(support_lasso)
print(a)
if(a<=p1){
  lam.min=lam.min/2
 }
}
}else{
  p0<-ncol(theta_snp)/2
  a<=0
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
fit<-FDR_selection_LM(data_list, bigdata=FALSE, lam.max=lam.max, lam.min=lam.min,nlambda=40,
                   alpha=0.1, tolerance=0.025, threshold=0.01,err=1e-5, delta=NULL)
Sys.time()-t1

support_FDR.T<-fit$support.T
b<-colnames(as.matrix(theta_snp[,support_FDR.T,drop=FALSE]))
cat("chr",k, "length", length(support_FDR.T), "support", b, "\n")

FDR_loss<-Loss_test(support_FDR.T, data_list)
FDR_loss
cat("lossT",(FDR_loss-test_loss)/test_loss,"\n")

support_FDR.F<-fit$support.F
b<-colnames(theta_snp[,support_FDR.F,drop=FALSE])
cat("chr",k, "length", length(support_FDR.F), "support", b, "\n")

FDR_loss<-Loss_test(support_FDR.F, data_list)
FDR_loss
cat("lossF",(FDR_loss-test_loss)/test_loss,"\n")
######## save result of selected theta_snps

theta_T<-cbind(theta_T,theta_snp[,support_FDR.T,drop=FALSE])

theta_F<-cbind(theta_F,theta_snp[,support_FDR.F,drop=FALSE])

rm(data_list,theta_snp,theta_train,theta_test,theta_val)
}

name<-paste0(geno_name,"_theta_T[1]",".RData")
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
FDR_loss<-Loss_test(support_FDR,data_list)
cat("loss_all_T",(FDR_loss-test_loss)/test_loss,"\n")


name<-paste0(geno_name,"_theta_F[1]",".RData")
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
FDR_loss<-Loss_test(support_FDR,data_list)
cat("loss_all_F",(FDR_loss-test_loss)/test_loss,"\n")

save.image("workspace[1].RData")
############BHq p_value###

# t1<-Sys.time()
# support_BHq<-BH_selection(data_list,family="binomial", alpha = 0.1, lambda=3e-3)$selected_snps
# support_BHq
# Sys.time()-t1
# 
# AUC_test(support_BHq, data_list, lasso=TRUE, lambda=3e-3, delta=1e-4)
# 


