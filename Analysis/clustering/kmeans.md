1. `slope1 = [p(Aa)-p(AA)]/(1-0)`
2. `slope2 = [p(aa)-p(Aa)]/(2-1)`
3. `ratio = slope2 / slope1` 
4. for SNPs with `ratio > 10, or ratio < -10`, put them into one group (might be recessive)
5. for remaining SNPs, run `kmeans`
