library(dplyr)

pheno <- "hypertension"

## functions for loading selected SNP

func_SNP_swap <- function(vec){
  # swap 0 and 2 in .raw
  vec[vec==0] <- 9
  vec[vec==2] <- 0
  vec[vec==9] <- 2
}

func_process_raw <- function(raw, SNP_list, SNP_complete){
  ## replace NA with 0
  raw <- as.matrix(raw)
  raw[is.na(raw)] <- 0
  colnames(raw) <- SNP_complete$SNP
  
  ## check mis-matched effective allele
  SNP_select <- inner_join(SNP_list, SNP_complete, by="SNP")
  SNP_select$match <- (SNP_select$effect_allele == SNP_select$effect_AoU)
  
  ## select SNPs from .raw
  raw_select <- raw[, SNP_select$SNP]
  
  if (nrow(SNP_select) > 0){
    for (i in ncol(raw_select)){
      # if effective allele mis-matched, swap 0 and 2
      if (!SNP_select$match[i]) {
        raw_select[,i] <- func_SNP_swap(raw_select[,i])
      }
    }
    if (nrow(SNP_select) == 1){
      ## if only one SNP selected
      ## force the vector to be column matrix
      raw_select <- matrix(raw_select, ncol=1)
      colnames(raw_select) <- SNP_select$SNP[1]
    }
    return(raw_select)
  } else {
    return(NULL)
  }
}

# Load raw and convert theta (for each chr)

SNP_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/selected_SNP/",pheno,"_SNP_list")
SNP_list <- read.table(SNP_dir, header=T)
#print(head(SNP_list))

betas_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/selected_SNP/",pheno,"_betas_mat")
betas_mat <- read.table(betas_dir, header=T)
#print(head(betas_mat))

for (i in 1:22){
  #print(paste("chr",i,"start..."))
  # read chr_i.raw
  raw_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/Version2/mainFiles/chr",i,"_small.raw")
  raw <- read.table(raw_dir, header=T)
  IID <- raw[,2]
  raw <- raw[,-c(1,2,3,4,5,6)]
  ## get the SNP list of .raw of the AoU
  SNP_complete <- data.frame(SNP = gsub("_.*","",colnames(raw)),
                             effect_AoU = gsub(".*_","",colnames(raw)))
  ## process .raw
  raw <- func_process_raw(raw, SNP_list, SNP_complete)
  
  ## if there are selected SNPs in chr_i, convert and save theta
  if (!is.null(raw)){
    ## convert theta chr_i
    df_selected <- data.frame(SNP=colnames(raw))
    betas_mat_temp <- inner_join(df_selected, betas_mat, by="SNP")
    theta <- convert_theta(raw, betas_mat_temp)
    
    ## save theta chr_i
    theta_dir <- paste0("/home/jupyter/workspaces/hypertensionukbb/AB_PRS_AoU/theta_SNP/",pheno,"_chr",i,"_theta")
    write.table(theta, theta_dir,quote=F,row.names=F,col.names=T)
  }
}




