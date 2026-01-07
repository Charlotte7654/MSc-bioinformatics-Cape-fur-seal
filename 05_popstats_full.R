# House keeping ####
#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts")

# Generate citations ####
#For R
citation("ade4")
citation("adegenet")

#For R Studio
if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")
library(rstudioapi)

versionInfo()$citation

#For packages
#citation("pkgname")

# Install and load packages ####
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)

BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

if (!require("hierfstat", quietly = TRUE))
  install.packages("hierfstat")
library(hierfstat)

if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)

if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

# Basic stats ####

# read vcf 
vcf_full <- read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_full.vcf")
genind_full <- vcfR2genind(vcf_full)

#used excel to add commas 
populations <- c(
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point"
)

length(populations) == nInd(genind_full) # should be true 

genind_full@pop <- as.factor(populations)

# convert genind obj to hierfstat object 
hf_full <- genind2hierfstat(genind_full)
is.genind(genind_full)
#TRUE

#Takes long
basic_stats_full <- basic.stats(hf_full)

# Save the object to a file on your computer
save(basic_stats_full, file = "basic_stats_full_data.RData")

# Load the saved object from the file
#load("basic_stats_full_data.RData")

#write a standard error function with NA omitted 
#se <- function(x) sd(x)/sqrt(length(x)) 
#produces same outputs in outlier file as funciton above 
se <- function(x) {
  x <- na.omit(x)            # remove NAs
  sd(x) / sqrt(length(x))    # now safe
}

#print overall stats 
basic_stats_full
#$overall
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.3137  0.2808  0.2808  0.0000  0.2808  0.0001  0.0002  0.0002 -0.1175  0.0001 

#Ho ####
#mean observed heterozygosities

Ho(hf_full) #takes long
#Cape Cross     False Bay      Kleinzee Lambert's Bay Pelican Point 
#    0.3143964     0.3134681     0.3150355     0.3103849     0.3153860 

Ho<-basic_stats_full[["Ho"]]
#View(Ho)

#Prep data for se calulation 
Ho_CC<-Ho[,1,drop=FALSE]
Ho_PP<-Ho[,5,drop=FALSE]
Ho_KZ<-Ho[,3,drop=FALSE]
Ho_LB<-Ho[,4,drop=FALSE]
Ho_FB<-Ho[,2,drop=FALSE]

#se calculation
se(Ho_CC) #0.0009445264
se(Ho_PP) #0.000944186
se(Ho_KZ) #0.0009446212
se(Ho_LB) #0.0009476427
se(Ho_FB) #0.0009708063

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Ho_CC, main="Histogram of Ho_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_CC)
qqline(Ho_CC, col = "red")
#CC not normal

hist(Ho_PP, main="Histogram of Ho_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_PP)
qqline(Ho_PP, col = "red")
#PP not normal

hist(Ho_KZ, main="Histogram of Ho_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_KZ)
qqline(Ho_KZ, col = "red")
#KZ not normal

hist(Ho_LB, main="Histogram of Ho_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_LB)
qqline(Ho_LB, col = "red")
#LB not normal

hist(Ho_FB, main="Histogram of Ho_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_FB)
qqline(Ho_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Ho_df <- as.data.frame(Ho)

Ho_long <- Ho_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Ho") # Name of the new 'Ho' column
# View the reshaped data
View(Ho_long)

kruskal.test(Ho ~ factor(population), data = Ho_long)
#	Kruskal-Wallis rank sum test
#data:  Ho by factor(population)
#Kruskal-Wallis chi-squared = 23.409, df = 4, p-value = 0.0001049

#Posthoc test for Kruskal-Wallis (Dunn's test)
# Perform Dunn's test
dunn_results_Ho <- dunn.test(Ho_long$Ho, Ho_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment

#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 23.4091, df = 4, p-value = 0

#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
#---------+--------------------------------------------
#False Ba |   1.962962
#         |     0.0414
#         |
#Kleinzee |  -0.547277  -2.510239
#         |     0.3245    0.0121*
#         |
#Lambert' |   3.107663   1.144701   3.654941
#|    0.0031*     0.1802    0.0006*
#  |
#  Pelican  |  -0.808689  -2.771651  -0.261411  -3.916352
#|     0.2617    0.0070*     0.3969    0.0004*
#  alpha = 0.05
#Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Ho_df <- data.frame(
  Comparison = dunn_results_Ho$comparisons,  # Pairwise comparisons
  Z = dunn_results_Ho$Z,                     # Z-statistics
  P.adjusted = dunn_results_Ho$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Ho_df)
#Comparison          Z   P.adjusted
#1         Cape Cross - False Bay  1.9629621 0.0413754820 # 
#2          Cape Cross - Kleinzee -0.5472777 0.3245488511
#3           False Bay - Kleinzee -2.5102398 0.0120649186 #
#4     Cape Cross - Lambert's Bay  3.1076638 0.0031428743 #
#5      False Bay - Lambert's Bay  1.1447017 0.1802376768 
#6       Kleinzee - Lambert's Bay  3.6549415 0.0006431015 #
#7     Cape Cross - Pelican Point -0.8086892 0.2616837138
#8      False Bay - Pelican Point -2.7716513 0.0069715934 #
#9       Kleinzee - Pelican Point -0.2614115 0.3968876117
#10 Lambert's Bay - Pelican Point -3.9163530 0.0004494930 #

# Significant differences:CC-FB, FB-KZ, CC-LB, KZ-LB, FB-PP, LB-PP 

#Hs ####
#mean gene diversities within population Hs
Hs(hf_full) #Takes long
#Cape Cross     False Bay      Kleinzee Lambert's Bay Pelican Point 
#    0.2804888     0.2811684     0.2815961     0.2790833     0.2814487 

Hs<-basic_stats_full[["Hs"]]

#Prep data for se calulation 
Hs_CC<-Hs[,1,drop=FALSE]
Hs_PP<-Hs[,5,drop=FALSE]
Hs_KZ<-Hs[,3,drop=FALSE]
Hs_LB<-Hs[,4,drop=FALSE]
Hs_FB<-Hs[,2,drop=FALSE]

se(Hs_CC) 
se(Hs_PP) 
se(Hs_KZ) 
se(Hs_LB) 
se(Hs_FB) 

> se(Hs_CC) 
[1] 0.0006440222
> se(Hs_PP) 
[1] 0.0006407662
> se(Hs_KZ) 
[1] 0.0006432898
> se(Hs_LB) 
[1] 0.000652455
> se(Hs_FB) 
[1] 0.0006639341

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Hs_CC, main="Histogram of Hs_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_CC)
qqline(Hs_CC, col = "red")
#CC not normal

hist(Hs_PP, main="Histogram of Hs_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_PP)
qqline(Hs_PP, col = "red")
#PP not normal

hist(Hs_KZ, main="Histogram of Hs_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_KZ)
qqline(Hs_KZ, col = "red")
#KZ almost normal

hist(Hs_LB, main="Histogram of Hs_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_LB)
qqline(Hs_LB, col = "red")
#LB not normal

hist(Hs_FB, main="Histogram of Hs_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_FB)
qqline(Hs_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Hs_df <- as.data.frame(Hs)

Hs_long <- Hs_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Hs") # Name of the new 'Hs' column

# View the reshaped data
View(Hs)
View(Hs_long)

kruskal.test(Hs ~ factor(population), data = Hs_long)
#	Kruskal-Wallis rank sum test
#data:  Hs by factor(population)
#Kruskal-Wallis chi-squared = 10.565, df = 4, p-value = 0.03191

#PostHoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#Dunn's test
dunn_results_Hs <- dunn.test(Hs_long$Hs, Hs_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment

#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 10.5653, df = 4, p-value = 0.03
#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
#False Ba |  -2.213477
#         |     0.0672
#         |
#Kleinzee |  -1.259481   0.954019
#         |     0.2079     0.2125
#         |
#Lambert' |   0.748237   2.961700   2.007719
#|     0.2524    0.0153*     0.0745
#|
#  Pelican  |  -1.022268   1.191228   0.237213  -1.770506
#|     0.2190     0.1946     0.4062     0.0958
#
#alpha = 0.05
#Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Hs_df <- data.frame(
  Comparison = dunn_results_Hs$comparisons,  # Pairwise comparisons
  Z = dunn_results_Hs$Z,                     # Z-statistics
  P.adjusted = dunn_results_Hs$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Hs_df)
Comparison          Z P.adjusted
#1         Cape Cross - False Bay -2.2134772 0.06716192 
#2          Cape Cross - Kleinzee -1.2594819 0.20785633
#3           False Bay - Kleinzee  0.9540192 0.21254620
#4     Cape Cross - Lambert's Bay  0.7482378 0.25239819
#5      False Bay - Lambert's Bay  2.9617007 0.01529725 #
#6       Kleinzee - Lambert's Bay  2.0077197 0.07445515
#7     Cape Cross - Pelican Point -1.0222684 0.21903850
#8      False Bay - Pelican Point  1.1912282 0.19463667
#9       Kleinzee - Pelican Point  0.2372135 0.40624557
#10 Lambert's Bay - Pelican Point -1.7705062 0.09580357

#Comparisons that are significant: FB-LB

#FIS ####
#Inbreeding Coeficcient 
#subset data to find FIS using basic.stats
#poppr for popsub
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

FIS_CC<-
  popsub(
    genind_full,
    sublist = "Cape Cross",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_PP<-
  popsub(
    genind_full,
    sublist = "Pelican Point",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_KZ<-
  popsub(
    genind_full,
    sublist = "Kleinzee",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_LB<-
  popsub(
    genind_full,
    sublist = "Lambert's Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_FB<-
  popsub(
    genind_full,
    sublist = "False Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

#Find FIS using basic.stats
basic.stats_CC <-basic.stats(FIS_CC,diploid=TRUE,digits=5)
basic.stats_CC
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31440  0.28049  0.28049  0.00000      NaN      NaN  0.00000      NaN -0.12088      NaN 

FIS_CC #-0.12088

basic.stats_PP<-basic.stats(FIS_PP,diploid=TRUE,digits=5)
basic.stats_PP
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31539  0.28145  0.28145  0.00000      NaN      NaN  0.00000      NaN -0.12058      NaN 

FIS_PP #-0.12058

basic.stats_KZ <-basic.stats(FIS_KZ,diploid=TRUE,digits=5)
basic.stats_KZ
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31504  0.28160  0.28160  0.00000      NaN      NaN  0.00000      NaN -0.11874      NaN 

FIS_KZ #-0.11874

basic.stats_LB <-basic.stats(FIS_LB,diploid=TRUE,digits=5)
basic.stats_LB
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31039  0.27909  0.27909  0.00000      NaN      NaN  0.00000      NaN -0.11216      NaN 

FIS_LB #-0.11216

basic.stats_FB <-basic.stats(FIS_FB,diploid=TRUE,digits=5)
basic.stats_FB

#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31347  0.28117  0.28117  0.00000      NaN      NaN  0.00000      NaN -0.11488      NaN 

FIS_FB #-0.11488

# I did this in the HPC 
#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=genind_full,nboot=1000,quant=c(0.025,0.975),diploid=TRUE,dig=4)
View(populations.structure_boot_fis)
populations.structure_boot_fis

View(populations.structure_boot_fis)
pops <- levels(genind_outlier@pop)
pops

??boot.ppfis hierfstat

#$call
#boot.ppfis(dat = genind_full, nboot = 1000, quant = c(0.025,
#                                                      0.975), diploid = TRUE, dig = 4)
#$fis.ci
#ll      hl
#1 -0.1244 -0.1169 CC
#2 -0.1188 -0.1110 FB
#3 -0.1221 -0.1152 KZ
#4 -0.1161 -0.1082 LB 
#5 -0.1243 -0.1169 PP 

#NB check order of pops ! 
# HPC done 

#Subset datasets to prepare to test nomality 
#Can probably also use basic.stats_PP$FIS
FIS<-basic_stats_full[["Fis"]]
View(FIS)
View(basic_stats_full)

#Prep data for se calulation 
FIS_CC<-FIS[,1,drop=FALSE]
FIS_PP<-FIS[,5,drop=FALSE]
FIS_KZ<-FIS[,3,drop=FALSE]
FIS_LB<-FIS[,4,drop=FALSE]
FIS_FB<-FIS[,2,drop=FALSE]

#se calculation
se(FIS_CC) #NA
se(FIS_PP) #NA
se(FIS_KZ) #NA
se(FIS_LB) #NA
se(FIS_FB) #NA

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(FIS_CC, main="Histogram of FIS_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_CC)
qqline(FIS_CC, col = "red")
#CC not normal

hist(FIS_PP, main="Histogram of FIS_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_PP)
qqline(FIS_PP, col = "red")
#PP not normal

hist(FIS_KZ, main="Histogram of FIS_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_KZ)
qqline(FIS_KZ, col = "red")
#KZ not normal

hist(FIS_LB, main="Histogram of FIS_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_LB)
qqline(FIS_LB, col = "red")
#LB not normal

hist(FIS_FB, main="Histogram of FIS_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_FB)
qqline(FIS_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
FIS_df <- as.data.frame(FIS)

FIS_long <- FIS_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "FIS") # Name of the new 'FIS' column
# View the reshaped data
View(FIS_long)

kruskal.test(FIS ~ factor(population), data = FIS_long)
#Kruskal-Wallis rank sum test
#data:  FIS by factor(population)
#Kruskal-Wallis chi-squared = 5.9031, df = 4, p-value = 0.2065

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#DON'T NEED POST-HOC (p>0.05)
# Perform Dunn's test
dunn_results_FIS <- dunn.test(FIS_long$FIS, FIS_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment


# Create a data frame from the results
dunn_results_FIS_df <- data.frame(
  Comparison = dunn_results_FIS$comparisons,  # Pairwise comparisons
  Z = dunn_results_FIS$Z,                     # Z-statistics
  P.adjusted = dunn_results_FIS$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_FIS_df)

#No comparisons are significant. No differences! 

#FST ####
#FST from Stamp 
if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

# Dependancies for dartR (not installing without them)
BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)


#full
populations.structure.full.gl<-gi2gl(genind_full)
save(populations.structure.full.gl, file = "populations.structure.full.gl.RData")

stamppFst_values <- stamppFst(populations.structure.gl, nboots = 100, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

# Save the object to a file on your computer
save(stamppFst_values, file = "stamppFst_values_full.RData")

# Load the saved object from the file
#load("stamppFst_values.RData")

options(max.print = 100)
stamppFst_values
View(stamppFst_values)

FST_pvalues<-stamppFst_values[["Pvalues"]]
FST_pvalues

#HDW (have not done) ####

# Read the VCF file
if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

#Check full to see if same as Populations
#data <- read.vcf("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.snps.vcf")
populations.snps.vcf <- read.vcf("populations.snps.vcf")
populations.snps.vcfR <- read.vcfR("populations.snps.vcf")

# Convert VCF to genotype matrix
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#.gen file was unavailable so I used the vcf
genind_full <- vcfR2genind(populations.snps.vcfR)

# Perform HWE test on the genind object using the `pegas` package
if (!require("pegas", quietly = TRUE))
  install.packages("pegas")
library(pegas)

HWE_results <- hw.test(genind_full)
# Write the data frame to a CSV file
write.csv(hwe_results, file = "HWE_results_full.csv", row.names = FALSE)

#Now full
#data <- read.vcf("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.snps.vcf")
populations.snps.netural.vcf <- read.vcf("populations.snps.full.vcf")

if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

populations.snps.netural.vcfR <- read.vcfR("populations.snps.full.vcf")

populations.snps.netural.vcfR <- read.vcfR("~/msc/00_info/populations.snps.full.vcf")

# Convert VCF to genotype matrix
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#.gen file was unavailable so I used the vcf
genind_full <- vcfR2genind(populations.snps.netural.vcfR)

# Perform HWE test on the genind object using the `pegas` package
if (!require("pegas", quietly = TRUE))
  install.packages("pegas")
library(pegas)

HWE_results_netural <- hw.test(genind_full)
# Write the data frame to a CSV file

write.csv(hwe_results_full, file = "HWE_results_full.csv", row.names = FALSE)
--> Run in hpc



