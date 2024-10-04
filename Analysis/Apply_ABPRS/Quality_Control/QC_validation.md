# QC of imputed PMBB and eMERGE data for AB-PRS validation

* Selected individuals who meet the following criteria 
  * Unrelated individual ($\hat{\pi} \geq 0.25$)
  * Recent European ancestry (ancestry-stratified)
  * Self-reported and inferred genetic sex match
  * Non-missing genotype or phenotype data
* Removed SNPs that meet the following critetia
  * MAF < 0.01
  * HWE p-value < 10
  * minimac $r^2$ < 0.3
  * Genotype missing rate > 0.1
  * PLINK `--hard-call-threshold 0.1`
