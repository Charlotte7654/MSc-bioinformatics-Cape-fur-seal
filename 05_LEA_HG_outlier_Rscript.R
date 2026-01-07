library(LEA)

snmf_outlier <- snmf("~/ch2/05_LEA/outlier/03.65_HG_outlier.geno",
K = 1:8, entropy = TRUE,
repetitions = 100,
project = "new",
alpha = 100)
