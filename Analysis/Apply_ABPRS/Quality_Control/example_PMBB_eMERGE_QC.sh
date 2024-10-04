plink2 \
    --pfile $input_prefix \
    --minimac3-r2-filter 0.3 2 \
    --hard-call-threshold 0.1 \
    --make-bed \
    --out ${output_prefix} \
    --threads ${plink_threads} \
    --memory ${plink_memory}

# perform basic QC
plink2 \
    --bfile $output_prefix \
    --maf 0.01 \
    --geno 0.1 \
    --hwe 1e-10 \
    --make-bed \
    --out ${output_prefix}_qc \
    --threads ${plink_threads} \
    --memory ${plink_memory}

# extract pruned SNPs and create a new dataset
plink2 \
    --bfile ${output_prefix}_qc \
    --extract range data/SNP_list.lifted \
    --remove /static/PMBB/PMBB-Release-2020-2.0/Imputed/Relatedness/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_relateds_droplist.txt \
    --recode A \
    --out ${output_prefix}_small \
    --threads ${plink_threads} \
    --memory ${plink_memory}
