library(LEA)

snmf_neutral <- snmf("~/ch2/05_LEA/neutral/03.65_HG_neutral.geno",
K = 1:8, entropy = TRUE,
repetitions = 100,
project = "new",
alpha = 100)
