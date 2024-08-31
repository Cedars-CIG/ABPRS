library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

print("loading icd10 data...")
icd10 <- fread(args[1],header=T,quote="")
icd10 <- as.data.frame(icd10)
print("loading icd10 data done.")


pheno_name <- c("Hypertension", "T2D", "AD", "Breast_Cancer")
code <- c("401.1", "250.2", "290.11", "174.1")

phenotypes <- icd10[, code]

print("get eid")
id <- icd10[,1]
phenotypes <- cbind(id, phenotypes)
colnames(phenotypes) <- c("eid", pheno_name)

print("write.csv...")
write.csv(phenotypes, file=args[2], quote=F, row.names=F, col.names=T)
print("write.csv done.")