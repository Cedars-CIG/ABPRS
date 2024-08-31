library(data.table)
library(dplyr)
library(PheWAS)

####################
##
## Discrete Phenotype
##
####################

icd10_prep_func <- function(dir_source, dir_icd10){
    ######### Response

    ### PHENO ###

    ## ICD 10 filter

    ## originl data
    print("loading source data")
    pheno<-fread(dir_source,header=T,quote="")

    print("getting primary icd10")
    ## 41202 primary diagnosis of icd10
    primary_icd10<-c()
    for (i in 0:78){
    name<-paste0("41202-0.",i)
    primary_icd10<-c(primary_icd10,name)
    }

    print("getting secondary icd10")
    ## 41204 secondary diagnosis of icd10
    secondary_icd10<-c()
    for (i in 0:187){
    name<-paste0("41204-0.",i)
    secondary_icd10<-c(secondary_icd10,name)
    }

    all_icd10 <- c(primary_icd10,secondary_icd10)

    ## select field 41202 and 41204 from pheno data
    icd10_pheno <- pheno[,all_icd10,with=F]

    ## create matrix contains id and pheno
    eid <- pheno[,"eid"]
    icd10_pheno <- cbind(eid,icd10_pheno)

    ## as data frame
    icd10_pheno <- as.data.frame(icd10_pheno)

    print("write table...")
    write.table(icd10_pheno, file=dir_icd10, quote=F, row.names=F, col.names=T, sep = ",")
    print("write table done.")

    #print(icd10_pheno[1:10, 1:10])
}

icd10_count_func <- function(dir_in, dir_out){
    print("loading original icd10 table")
    icd10<-fread(dir_in,header=T,sep=",")
    icd10 <- as.data.frame(icd10)
    #print(head(icd10))

    # get id
    id <- icd10$eid
    icd10 <- icd10[,-c(1)]

    # final data frame col order
    # match with the example
    col_order <- c("id", "vocabulary_id", "code", "count")

    # 1st individual
    vec <- as.character(icd10[1,])
    vec <- vec[vec != ""]
    df <- as.data.frame(table(vec))
    colnames(df) <- c("code", "count")
    df$id <- rep(id[1], nrow(df))
    df$vocabulary_id <- rep("ICD10", nrow(df))
    df <- df[, col_order]

    print("start getting final table")
    print(paste("total number of individuals:", nrow(icd10)))
    # remaining individuals
    for (i in 2:nrow(icd10)){

        if (i %% 5000 == 0){print(paste("individual", i, "done."))}
        vec <- as.character(icd10[i,])
        vec <- vec[vec != ""]
        if (length(vec)>0) {
            df_df_covariates <- as.data.frame(table(vec))
            colnames(df_df_covariates) <- c("code", "count")
            df_df_covariates$id <- rep(id[i], nrow(df_df_covariates))
            df_df_covariates$vocabulary_id <- rep("ICD10", nrow(df_df_covariates))
            df_df_covariates <- df_df_covariates[, col_order]
            df <- rbind(df, df_df_covariates)
        }
    }
    print("write.table...")
    write.table(df, file=dir_out, quote=F, row.names=F, col.names=T)
    print("write.table done.")
}


icd10_code_modify <- function(code){
    part1 <- substr(code, 1, 3)
    part2 <- substr(code, 4, 4)
    code_out <- paste0(part1, ".", part2)
    return(code_out)
}

icd10_phenotype_func <- function(dir_in, dir_out){
    print("loading icd10 count")
    icd10_count <- fread(dir_in)
    icd10_count <- as.data.frame(icd10_count)
    icd10_count$vocabulary_id <- factor(icd10_count$vocabulary_id)
    #icd10_count$vocabulary_id <- factor(rep("ICD10CM", nrow(icd10_count)))
    print("total number of rows")
    print(nrow(icd10_count))
    print("total number of NA in $count")
    print(sum(is.na(icd10_count$count)))
    print("data type")
    print(str(icd10_count))

    print("create phenotypes...")
    phenotypes <- createPhenotypes(icd10_count, 
                                   min.code.count = 1,
                                   vocabulary.map = PheWAS::phecode_map_icd10)
    phenotypes <- as.data.frame(phenotypes)
    print("create phenotypes done.")

    print("original phenotype with TRUE/FALSE")
    print(phenotypes[1:5, 1:5])

    print("running lapply")
    phenotypes[] <- lapply(phenotypes, as.integer)
    print("new phenotype with 0/1")
    print(phenotypes[1:5, 1:5])

    print("write.csv...")
    write.csv(phenotypes, file=dir_out, quote=F, row.names=F, col.names=T)
    print("write.csv done.")
}







