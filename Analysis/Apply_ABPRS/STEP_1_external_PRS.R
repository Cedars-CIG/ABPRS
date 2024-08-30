library(dplyr)

pheno <- "hypertension"

## load GWAS (copied from cedars server)
GWAS_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/external_GWAS/",pheno,"_external_GWAS")
GWAS <- read.table(GWAS_dir, header=T)
#print(head(GWAS))

## PRS calculation

func_SNP_swap <- function(vec){
  # swap 0 and 2 in .raw
  vec[vec==0] <- 9
  vec[vec==2] <- 0
  vec[vec==9] <- 2
}

func_process_raw <- function(raw, GWAS, SNP_complete){
  ## replace NA with 0
  raw <- as.matrix(raw)
  raw[is.na(raw)] <- 0
  colnames(raw) <- SNP_complete$SNP
  
  ## check mis-matched effective allele
  SNP_select <- inner_join(GWAS, SNP_complete, by="SNP")
  SNP_select$match <- (SNP_select$effect_allele == SNP_select$effect_AoU)
  
  ## select SNPs from .raw
  raw_select <- raw[, SNP_select$SNP]
  
  if (ncol(raw_select) > 0){
    for (i in ncol(raw_select)){
      # if effective allele mis-matched, swap 0 and 2
      if (!SNP_select$match[i]) {
        raw_select[,i] <- func_SNP_swap(raw_select[,i])
      }
    }
    return(list(raw_select, SNP_select))
  } else {
    return(list(raw_select, SNP_select))
  }
}

func_prs <- function(rslt){
  raw <- rslt[[1]]
  gwas <- rslt[[2]]
  prs <- raw %*% gwas$beta
  return(prs)
}


#print(paste("chr",1,"start..."))
# read chr_1.raw
raw_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/Version2/mainFiles/chr",1,"_small.raw")
raw <- read.table(raw_dir, header=T)
IID <- raw[,2]
raw <- raw[,-c(1,2,3,4,5,6)]
## get the SNP list of .raw of the AoU
SNP_complete <- data.frame(SNP = gsub("_.*","",colnames(raw)),
                           effect_AoU = gsub(".*_","",colnames(raw)))
rslt <- func_process_raw(raw, GWAS, SNP_complete)
prs <- func_prs(rslt)

for (i in 2:22){
  #print(paste("chr",i,"start..."))
  # read chr_i.raw
  raw_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/Version2/mainFiles/chr",i,"_small.raw")
  raw <- read.table(raw_dir, header=T)
  raw <- raw[,-c(1,2,3,4,5,6)]
  ## get the SNP list of .raw of the AoU
  SNP_complete <- data.frame(SNP = gsub("_.*","",colnames(raw)),
                             effect_AoU = gsub(".*_","",colnames(raw)))
  rslt <- func_process_raw(raw, GWAS, SNP_complete)
  
  ## if there are selected SNPs in chr_i
  if (nrow(rslt[[2]]) > 0){
    prs_temp <- func_prs(rslt)
    prs <- prs + prs_temp
  } 
}


## save prs dataframe
df_rslt <- data.frame(IID=IID, prs=prs)
#print(dim(df_rslt))
#print(head(df_rslt))

prs_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/external_PRS/",pheno,"_external_PRS.csv")
write.csv(df_rslt, prs_dir, row.names = FALSE)

