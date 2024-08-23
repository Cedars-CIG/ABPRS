Clustering algorithm:
1. `slope1 = [p(Aa)-p(AA)]/(1-0)`
2. `slope2 = [p(aa)-p(Aa)]/(2-1)`
3. `ratio = slope2 / slope1` 
4. For SNPs with `ratio > 10, or ratio < -10`, put them into one group (might be recessive)
5. For remaining SNPs, run `kmeans`

Several details:
1. Based on the sum of squared plots, chose 7 as number of cluster for `kmeans`
2. So in total, we have 7+1 clusters
3. For discrete phenotype, some of the SNPs missing SNP names
