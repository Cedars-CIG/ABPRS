- For participants
    1. Filter white British
    2. Match genetic and report sex
    3. Check heter < mu+3*sd, heter > mu-3*sd
    4. Remove NA Phenotype (Hypertension / T2D / HDL / LDL…)
    5. Drop related participants
    6. Remove participants that DO NOT have genotype data
- For SNPs
    1. Get SNP list based on Neale GWAS (for --extract flag)
        1. Only select SNPs with p-val < 0.05
        2. Remove ambiguous SNPs (using R)
    2. Get A1 list from Neale GWAS
    3. Remove duplicates using PLINK: --exclude chr_X.dup
        1. Do NOT remove duplicates when generate .bed and .bim
        2. Get .dup based on .bim file
    4. Add hwe flag in converting .raw
