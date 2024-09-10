Datasets
- Phenotype: please refer to the phenotype generation process in `Phenotype` folder;
- Genotype:
    - `N by P` matrix of 0/1/2 coding. We used `.raw` from `PLINK`, but you can also use other ways to generate this `N by P` matrix.
    - `N` is the number of individuals, and `P` is the number of SNPs.
    - For external PRS calculation: all GWAS information (including `chr`, `bp`, `rs_xxx` `effect_allele`, `effect_size`) is listed in `eternal_GWAS` folder.
    - For AB-PRS analysis: the information (including `chr`, `bp`, `rs_xxx` `effect_allele`) of seletected SNP is listed in `selected_SNP` folder.

Analysis

Please follow STEP 1 to STEP 3, and refer to `.ipynb` for more details
