library(LEA)

snmf_full <- snmf("~/ch2/05_LEA/full/03.65_HG_full.geno",
K = 1:8, entropy = TRUE,
repetitions = 100,
project = "new",
alpha = 100)
