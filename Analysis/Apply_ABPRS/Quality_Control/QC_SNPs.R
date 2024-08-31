library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

print(args[1])

print("loading data...")
GWAS <- read.table(args[1], header=T)
print("loading data done.")

print("# SNPs, total")
print(nrow(GWAS))

##### pval < 0.05
pval <- as.numeric(args[2])
GWAS <- GWAS[(GWAS$pval < pval), ]
print(paste("# SNPs, keep pval <", pval))
print(nrow(GWAS))

##### remove ambiguous
GWAS <- GWAS[!(GWAS$allele_pair %in% c("AT", "TA", "CG", "GC")), ]
print("# SNPs, after remove ambiguous")
print(nrow(GWAS))

SNP <- GWAS %>% select("SNP")
A1 <- GWAS %>% select("SNP", "minor_allele")

write.table(SNP, file=args[3], quote=F, row.names=F, col.names=T)
write.table(A1, file=args[4], quote=F, row.names=F, col.names=T)
print("write.csv done.")