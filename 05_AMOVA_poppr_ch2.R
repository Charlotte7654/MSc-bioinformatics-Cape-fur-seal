#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts")

#AMOVA ####
# Load required libraries if not already loaded
if (!require("vcfR", quietly = TRUE)) install.packages("vcfR")
library(vcfR)
if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
library(adegenet)
if (!require("poppr", quietly = TRUE)) install.packages("poppr")
library(poppr)
if (!require("pegas", quietly = TRUE)) install.packages("pegas")
library(pegas)

# think I need ade4 package too! 

# proportion of missing data per locus
miss_prop <- function(full.genind) {
  n_loci <- nLoc(full.genind)
  sapply(1:n_loci, function(i) {
    # mean proportion of missing genotypes for locus i
    mean(is.na(full.genind@tab[, i]))
  })
}

missing_per_locus <- miss_prop(full.genind)
n_retained <- sum(missing_per_locus <= 0.05)
n_retained
20% 20742
15% 11864
10% 5194
7% 1652
6% 934
5% 426

# proportion of missing data per locus
miss_prop_neutral <- function(neutral.genind) {
  n_loci_neutral <- nLoc(neutral.genind)
  sapply(1:n_loci_neutral, function(i) {
    # mean proportion of missing genotypes for locus i
    mean(is.na(neutral.genind@tab[, i]))
  })
}

missing_per_locus_neutral <- miss_prop(neutral.genind)
n_retained_neutral <- sum(missing_per_locus_neutral <= 0.05)
n_retained_neutral
20% 20664
15% 11842
10% 5190
7% 1648
6% 932
5% 426

missing_per_locus_outlier <- miss_prop(outlier.genind)
n_retained_outlier <- sum(missing_per_locus_outlier <= 0.31)
n_retained_outlier
30% 166
20% 84
15% 26
10% 6
7% 4
6% 2
5% 0 

#from AMOVA: 158 loci contained missing values greater than 30% 
num = 324-158
num


# FULL #### 
# Load genind, if it exists####
load("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/full.genind.RData")

# SKIP (load VCF and create genind) ####
full.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_full.vcf")
full.genind <- vcfR2genind(full.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#individual	population
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Assign population factor to genind "pop" based on metadata (full name)
pop(full.genind)<- as.factor(popmap$population)
pop(full.genind)

# Assign country and population factors to genind strata based on metadata (full name)
#strata(full.genind) <- data.frame(country = popmap$country, population = popmap$population)

# Assign  population factors to genind strata based on metadata (full name)
strata(full.genind) <- data.frame(population = popmap$population)

# Run AMOVA using the poppr package ####
# no need to use nested model. 
#SKIP First between countries, then between pops within countries, then between pops. 
#c.p_poppr_amova_full <- poppr.amova(full.genind, ~ country / population)
#print(c.p_poppr_amova_full)

# SKIP (frist country, then population)
#cp_poppr_amova_full <- poppr.amova(full.genind, ~ country + population)
#print(cp_poppr_amova_full)

# I just want to group by population (takes long!)
#c.p_poppr_amova_full <- poppr.amova(full.genind, ~ population)
#print(c.p_poppr_amova_full)

# NOT FILTERED 
c.p_poppr_amova_full <- poppr.amova(full.genind, ~ population)
print(c.p_poppr_amova_full)

c.p_poppr_amova_full_rand <- randtest(c.p_poppr_amova_full, nrepet = 100000)
c.p_poppr_amova_full_rand

?? randtest 
#from ade4

# NOT FILTERED (5% default missingness tolerance)
#$call
#ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

#$results
#Df    Sum Sq  Mean Sq
#Between samples  4  240.3516 60.08791
#Within samples  86 4735.2507 55.06105
#Total           90 4975.6023 55.28447

#$componentsofcovariance
#Sigma           %
#Variations  Between samples  0.2764846   0.4996331
#Variations  Within samples  55.0610543  99.5003669
#Total variations            55.3375390 100.0000000

#$statphi
#Phi
#Phi-samples-total 0.004996331

#Monte-Carlo test
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 0.2764846 
#Based on 100000 replicates
#Simulated p-value: 0.01716983 
#Alternative hypothesis: greater 
#Std.Obs  Expectation     Variance 
#2.1614351257 0.0002163359 0.0163372044 

# FILTERED 
c.p_poppr_amova_full_filtered <- poppr.amova(full.genind, 
                                             ~ population, 
                                             missing = "loci", 
                                             cutoff = 0.06)
print(c.p_poppr_amova_full_filtered)

c.p_poppr_amova_full_filtered_rand <- randtest(c.p_poppr_amova_full_filtered, nrepet = 100000)
c.p_poppr_amova_full_filtered_rand

# FILTERED 
#$call
#ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

#$results
#Df     Sum Sq  Mean Sq
#Between samples  4   516.7616 129.1904
#Within samples  86 10520.4424 122.3307
#Total           90 11037.2040 122.6356

#$componentsofcovariance
#Sigma           %
#Variations  Between samples   0.377292   0.3074714
#Variations  Within samples  122.330726  99.6925286
#Total variations            122.708018 100.0000000

#$statphi
#Phi
#Phi-samples-total 0.003074714

#Monte-Carlo test
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 0.377292 
#Based on 100000 replicates
#Simulated p-value: 0.03244968 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#1.881514701 0.001203752 0.039954377 

plot(c.p_poppr_amova_full_rand)

# Save AMOVA results to RDS file
saveRDS(c.p_poppr_amova_full, "/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts/poppr_amova_result_full.rds")

# Reload AMOVA results from RDS file
#poppr_amova_result_loaded <- readRDS("C:/Users/Charlotte/Documents/poppr_amova_result.rds")


# NEUTRAL #### 
# Load genind, if it exists####
load("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/neutral.genind.RData")

# SKIP (load VCF and create genind) ####
neutral.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_neutral.vcf")
neutral.genind <- vcfR2genind(neutral.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#individual	population
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Assign population factor to genind "pop" based on metadata (neutral name)
pop(neutral.genind)<- as.factor(popmap$population)
pop(neutral.genind)

# Assign country and population factors to genind strata based on metadata (neutral name)
#strata(neutral.genind) <- data.frame(country = popmap$country, population = popmap$population)

# Assign  population factors to genind strata based on metadata (neutral name)
strata(neutral.genind) <- data.frame(population = popmap$population)

# Run AMOVA using the poppr package ####
# no need to use nested model. 
#SKIP First between countries, then between pops within countries, then between pops. 
#c.p_poppr_amova_neutral <- poppr.amova(neutral.genind, ~ country / population)
#print(c.p_poppr_amova_neutral)

# SKIP (frist country, then population)
#cp_poppr_amova_neutral <- poppr.amova(neutral.genind, ~ country + population)
#print(cp_poppr_amova_neutral)

# I just want to group by population (takes long!)
#c.p_poppr_amova_neutral <- poppr.amova(neutral.genind, ~ population)
#print(c.p_poppr_amova_neutral)

c.p_poppr_amova_neutral <- poppr.amova(neutral.genind)
print(c.p_poppr_amova_neutral)

c.p_poppr_amova_neutral_rand <- randtest(c.p_poppr_amova_neutral, nrepet = 100000)
c.p_poppr_amova_neutral_rand

#NON FILTERED
#$call
#ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

#$results
#Df    Sum Sq  Mean Sq
#Between samples  4  240.3516 60.08791
#Within samples  86 4735.2507 55.06105
#Total           90 4975.6023 55.28447

#$componentsofcovariance
#Sigma           %
#Variations  Between samples  0.2764846   0.4996331
#Variations  Within samples  55.0610543  99.5003669
#Total variations            55.3375390 100.0000000

#Monte-Carlo test
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 0.2764846 
#Based on 100000 replicates
#Simulated p-value: 0.01697983 
#Alternative hypothesis: greater 
#Std.Obs   Expectation      Variance 
#2.1542951542 -0.0004009468  0.0165192497 

#FILTERED 
c.p_poppr_amova_neutral_filtered <- poppr.amova(neutral.genind, 
                                                ~ population, 
                                                missing = "loci", 
                                                cutoff = 0.06)
print(c.p_poppr_amova_neutral_filtered)

c.p_poppr_amova_neutral_filtered_rand <- randtest(c.p_poppr_amova_neutral_filtered, nrepet = 100000)
c.p_poppr_amova_neutral_filtered_rand

#FILTERED 
#$call
#ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

#$results
#Df     Sum Sq  Mean Sq
#Between samples  4   515.0982 128.7746
#Within samples  86 10490.2992 121.9802
#Total           90 11005.3974 122.2822

#$componentsofcovariance
#Sigma           %
#Variations  Between samples   0.3736985   0.3054242
#Variations  Within samples  121.9802228  99.6945758
#Total variations            122.3539214 100.0000000

#$statphi
#Phi
#Phi-samples-total 0.003054242
#Monte-Carlo test
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 0.3736985 
#Based on 100000 replicates
#Simulated p-value: 0.03264967 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#1.876361987 0.000592106 0.039539596 

plot(c.p_poppr_amova_neutral_rand)

# Save AMOVA results to RDS file
saveRDS(c.p_poppr_amova_neutral, "/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts/poppr_amova_result_neutral.rds")

# Reload AMOVA results from RDS file
#poppr_amova_result_loaded <- readRDS("C:/Users/Charlotte/Documents/poppr_amova_result.rds")

# OUTLIER #### 
# Load genind, if it exists####
load("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/outlier.genind.RData")

# SKIP (load VCF and create genind) ####
outlier.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_outlier.vcf")
outlier.genind <- vcfR2genind(outlier.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#individual	population
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Assign population factor to genind "pop" based on metadata (outlier name)
pop(outlier.genind)<- as.factor(popmap$population)
pop(outlier.genind)

# Assign country and population factors to genind strata based on metadata (outlier name)
#strata(outlier.genind) <- data.frame(country = popmap$country, population = popmap$population)

# Assign  population factors to genind strata based on metadata (outlier name)
strata(outlier.genind) <- data.frame(population = popmap$population)

# Run AMOVA using the poppr package ####
# no need to use nested model. 
#SKIP First between countries, then between pops within countries, then between pops. 
#c.p_poppr_amova_outlier <- poppr.amova(outlier.genind, ~ country / population)
#print(c.p_poppr_amova_outlier)

# SKIP (frist country, then population)
#cp_poppr_amova_outlier <- poppr.amova(outlier.genind, ~ country + population)
#print(cp_poppr_amova_outlier)

# I just want to group by population (takes long!)
#c.p_poppr_amova_outlier <- poppr.amova(outlier.genind, ~ population, missing = "ignore")
#print(c.p_poppr_amova_outlier)


# FILTERED 

c.p_poppr_amova_outlier_filtered <- poppr.amova(outlier.genind, 
                                                 ~ population, 
                                                 missing = "loci", 
                                                 cutoff = 0.30)
print(c.p_poppr_amova_outlier_filtered)

c.p_poppr_amova_outlier_filtered_rand <- randtest(c.p_poppr_amova_outlier_filtered, nrepet = 100000)
c.p_poppr_amova_outlier_filtered_rand


c.p_poppr_amova_outlier_filtered2 <- poppr.amova(outlier.genind, 
                                            ~ population, 
                                            missing = "loci", 
                                            cutoff = 0.50)
print(c.p_poppr_amova_outlier_filtered2)

c.p_poppr_amova_outlier_filtered_rand2 <- randtest(c.p_poppr_amova_outlier_filtered, nrepet = 100000)
c.p_poppr_amova_outlier_filtered_rand2


# FILTERED 
#$call
#ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

#$results
#Df     Sum Sq  Mean Sq
#Between samples  4   84.65981 21.16495
#Within samples  86 1129.82470 13.13750
#Total           90 1214.48451 13.49427

#$componentsofcovariance
#Sigma          %
#Variations  Between samples  0.4415222   3.251503
#Variations  Within samples  13.1374965  96.748497
#Total variations            13.5790187 100.000000

#Monte-Carlo test
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 0.4415222 
#Based on 100000 replicates
#Simulated p-value: 0.0003699963 
#Alternative hypothesis: greater 
#Std.Obs  Expectation     Variance 
#4.6802351314 0.0002829544 0.0088881830 


#plot(c.p_poppr_amova_outlier_rand)
plot(c.p_poppr_amova_outlier_filtered_rand)

# Save AMOVA results to RDS file
saveRDS(c.p_poppr_amova_outlier, "/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts/poppr_amova_result_outlier.rds")

# Reload AMOVA results from RDS file
#poppr_amova_result_loaded <- readRDS("C:/Users/Charlotte/Documents/poppr_amova_result.rds")

