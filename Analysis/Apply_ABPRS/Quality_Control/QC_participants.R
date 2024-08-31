library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

print(args[1])

if (args[1] %in% c("Hypertension", "T2D", "Breast_Cancer")){
    df_pheno <- read.csv("/common/lir2lab/Wanling/DATA/PheWAS/discrete/phenotype_discrete.csv", header=T)
}
if (args[1] == "AD_DNAnexus"){
    df_pheno <- read.csv("/common/lir2lab/Wanling/DATA/DNAnexus/phenotype/AD_DNAnexus_phenotype.csv", header=T)
}
if (args[1] %in% c("LDL", "HDL", "Cholesterol", "BMI")){
    df_pheno <- read.csv("/common/lir2lab/Wanling/DATA/UKBB/Phenotype/continuous/phenotype_continuous.csv", header=T)
}

df_cov <- read.csv(args[2], header=T)
df_cov <- as.data.frame(df_cov)
print("total number of participants")
print(nrow(df_cov))


df_IID_complete <- data.frame(IID = df_cov$eid)
df_IID <- data.frame(IID = df_cov$eid)

############################
## QC, drop related
############################

print('QC, drop related')

drop_related <- read.table(args[3],header=T)

## filter the IID, drop related
df_IID <- df_IID[which(!(df_IID$IID %in% drop_related$ID1)),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after drop related")
print(nrow(df_IID))

###############################
## QC Ancestry: white British
###############################

print('QC Ancestry: white British')
ancestry <- df_cov %>% select(all_of(paste0("Ancestry_",seq(1:3))))
ancestry <- as.data.frame(ancestry)
## return an true/false matrix/df of 1001(white british)
position <- apply(ancestry, 2, function(x) grepl("1001",x))

## true -> 1, na/false -> 0
ancestry[position] <- 1
ancestry[ancestry!=1] <- 0
ancestry[is.na(ancestry)] <- 0

## select all white british
ancestry <- as.data.frame(lapply(ancestry, as.numeric))
ancestry_status <- as.numeric(apply(ancestry,1,sum)>0)
ancestry_british <- cbind(df_IID_complete, ancestry_status)[which(ancestry_status==1),1]

## filter the IID of selected white british
df_IID <- df_IID[which(df_IID$IID %in% ancestry_british),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after white British")
print(nrow(df_IID))

###############################
## QC heter +/- 3sd
###############################

print("QC heter +/- 3sd")

heter <- df_cov %>% select("Heterozygosity")
heter <- as.data.frame(heter)

mu <- mean(heter$Heterozygosity, na.rm = TRUE)
sd <- sd(heter$Heterozygosity, na.rm = TRUE)

upper <- mu + 3 * sd
lower <- mu - 3 * sd

## return an true/false matrix/df within 3sd)
position <- apply(heter, 2, function(x) (lower <= x) & (x <= upper))

## true -> 1, na/false -> 0
heter[position] <- 1
heter[heter!=1] <- 0
heter[is.na(heter)] <- 0
heter <- as.data.frame(lapply(heter, as.numeric))
## get the valid id
heter_w3sd <- cbind(df_IID_complete, heter)[which(heter==1),1]

## filter the IID of selected id
df_IID <- df_IID[which(df_IID$IID %in% heter_w3sd),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after heter +/- 3sd")
print(nrow(df_IID))

###############################
## QC match sex
###############################

print("QC match sex")

sex <- df_cov %>% select(all_of(c("eid", "Sex_Genetic", "Sex_Report")))

sex <- sex[sex$Sex_Genetic == sex$Sex_Report, c(1,2)]
colnames(sex) <- c("IID", "Sex")

## filter the IID of matched sex
df_IID <- df_IID[which(df_IID$IID %in% sex$IID),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after match sex")
print(nrow(df_IID))

########################
## remove NA phenotype
########################

print("remove NA phenotype")
pheno_name <- args[1]
pheno <- df_pheno %>% select(all_of(c("eid", pheno_name)))
colnames(pheno) <- c("IID", "phenotype")

pheno <- pheno[!is.na(pheno$phenotype),]

## filter the IID of matched sex
df_IID <- df_IID[which(df_IID$IID %in% pheno$IID),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after remove NA phenotype")
print(nrow(df_IID))


########################
## inner join with .fam
########################

print("inner join with .fam")
fam <- read.table(args[4], header=F)
colnames(fam) <- c("IID", "V2", "V3", "V4", "V5", "V6")

## filter the IID of matched sex
df_IID <- df_IID[which(df_IID$IID %in% fam$IID),]
df_IID <- data.frame(IID = df_IID)

print("# participants, after inner join with .fam")
print(nrow(df_IID))


#####################################
## remove male (only breast cancer)
#####################################

if (pheno_name == "Breast_Cancer"){
    print("breast cancer, remove male")
    rm_male <- inner_join(df_IID, sex, by="IID")
    rm_male <- rm_male[rm_male$Sex == 0, ]
    df_IID <- data.frame(IID = rm_male$IID)

    print("# participants, after remove male (only breast cancer)")
    print(nrow(df_IID))
}


#################################
## Output Phenotype & Covariates
#################################

print("get phenotype, inner join with IID")
pheno <- inner_join(df_IID, pheno, by="IID")
print("# participants, phenotype")
print(nrow(pheno))

print("get covariates, inner join with IID")
cov_list <- c("eid", paste0("PC",seq(1:20)), "Age", "Batch")
cov <- df_cov %>% select(all_of(cov_list))
colnames(cov) <- c("IID",paste0("PC",seq(1:20)), "Age", "Batch")
cov <- inner_join(cov, sex, by="IID")

cov <- inner_join(df_IID, cov, by="IID")

print("# participants, covariates")
print(nrow(cov))

fam <- c(df_IID[,c(1,1)],0,0,0,0)

print("write.table start...")
write.table(df_IID,args[5],col.names=F,row.names=F,quote=F)
write.table(fam,args[6],col.names=F,row.names=F,quote=F)
write.table(pheno,args[7],col.names=T,row.names=F,quote=F)
write.table(cov,args[8],col.names=T,row.names=F,quote=F)
print("wrtie.table done.")

