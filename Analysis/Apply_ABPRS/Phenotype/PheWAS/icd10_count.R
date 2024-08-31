source("/common/lir2lab/Wanling/DATA/PheWAS/Rscript/PheWAS_library.R")

args <- commandArgs(TRUE)

icd10_count_func(args[1], args[2])
