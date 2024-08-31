library(data.table)
library(dplyr)

####################
##
## Continuous Phenotype
##
####################

pheno_continuous_func <- function(pheno, code){
    df <- pheno[, paste0(code, "-0.0"), with=F]
    df <- as.data.frame(df)
    return(df[,1])
}

args <- commandArgs(TRUE)

print("loading source data...")
pheno<-fread(args[1],header=T,quote="")
print("loading source data done.")


pheno_name <- c("LDL", "HDL", "Cholesterol", "BMI")
code <- c(30780, 30760, 30690, 21001)
#inst <- c(2, 2, 2, 4)

vec_list <- vector("list", 4)
for (i in 1:4){
    print(pheno_name[i])
    vec_list[[i]] <- pheno_continuous_func(pheno, code[i])
}

phenotypes <- as.data.frame(do.call(cbind, vec_list))  

print("get IID")
id <- pheno[,"eid"]
phenotypes <- cbind(id, phenotypes)
colnames(phenotypes) <- c("eid", pheno_name)

print("write.csv...")
write.csv(phenotypes, file=args[2], quote=F, row.names=F, col.names=T)
print("write.csv done.")