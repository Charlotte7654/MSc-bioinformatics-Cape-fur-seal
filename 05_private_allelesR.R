#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts")

if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

# FULL ####
# load VCF and create genind #
full.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_full.vcf")
full.genind <- vcfR2genind(full.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#sample	pop
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Ensure the columns used for matching are correctly named
names(popmap) <- c("individual", "population") # Standardize names if they differ from the file layout

# Assign population factor to genind "pop" based on metadata (full name)
pop(full.genind)<- as.factor(popmap$population)
pop(full.genind)

## levels(full.genind$pop) <- c("K", "PP", "CC", "M", "S", "LJP", "TS", "SR") #Based on order in popmap_ref_92

poppr_geneclone <- as.genclone(full.genind)

# Finding private alleles by population
private_alleles_result <- poppr::private_alleles(poppr_geneclone, level = "population")

# Convert the result to a data frame
private_alleles_df <- as.data.frame(private_alleles_result)
View(private_alleles_df)

#Don't run (print)
# Function to display private alleles by population
display_private_alleles <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    cat("\nPrivate alleles for population:", pop, "\n")
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print the private alleles and their count
    if (length(private_alleles_pop) > 0) {
      print(private_alleles_pop)
      cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
    } else {
      cat("No private alleles for population", pop, "\n")
    }
  }
}

#Don't run (print)
# Assuming private_alleles_result is your input data (e.g., allele frequencies or counts)
# Example of how to call the function
display_private_alleles(private_alleles_result)

#Run (counts)
# Function to display only the count of private alleles by population
display_private_allele_counts <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print only the count of private alleles
    cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
  }
}

#Run (counts)
# Example of how to call the function
display_private_allele_counts(private_alleles_result)

#Number of private alleles for population Cape Cross : 75 
#Number of private alleles for population False Bay : 58 
#Number of private alleles for population Kleinzee : 93 
#Number of private alleles for population Lambert's Bay : 72 
#Number of private alleles for population Pelican Point : 70 

# NEUTRAL #### 
# load VCF and create genind #
neutral.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_neutral.vcf")
neutral.genind <- vcfR2genind(neutral.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#sample	pop
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Ensure the columns used for matching are correctly named
names(popmap) <- c("individual", "population") # Standardize names if they differ from the file layout

# Assign population factor to genind "pop" based on metadata (neutral name)
pop(neutral.genind)<- as.factor(popmap$population)
pop(neutral.genind)

## levels(neutral.genind$pop) <- c("K", "PP", "CC", "M", "S", "LJP", "TS", "SR") #Based on order in popmap_ref_92

poppr_geneclone <- as.genclone(neutral.genind)

# Finding private alleles by population
private_alleles_result <- poppr::private_alleles(poppr_geneclone, level = "population")

# Convert the result to a data frame
private_alleles_df <- as.data.frame(private_alleles_result)
View(private_alleles_df)

#Don't run (print)
# Function to display private alleles by population
display_private_alleles <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    cat("\nPrivate alleles for population:", pop, "\n")
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print the private alleles and their count
    if (length(private_alleles_pop) > 0) {
      print(private_alleles_pop)
      cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
    } else {
      cat("No private alleles for population", pop, "\n")
    }
  }
}

#Don't run (print)
# Assuming private_alleles_result is your input data (e.g., allele frequencies or counts)
# Example of how to call the function
display_private_alleles(private_alleles_result)

#Run (counts)
# Function to display only the count of private alleles by population
display_private_allele_counts <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print only the count of private alleles
    cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
  }
}

#Run (counts)
# Example of how to call the function
display_private_allele_counts(private_alleles_result)
#Number of private alleles for population Cape Cross : 74 
#Number of private alleles for population False Bay : 57 
#Number of private alleles for population Kleinzee : 92 
#Number of private alleles for population Lambert's Bay : 71 
#Number of private alleles for population Pelican Point : 70 

# OUTLIER ####
# load VCF and create genind #
outlier.vcf <-read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_outlier.vcf")
outlier.genind <- vcfR2genind(outlier.vcf)

# Load metadata file containing individual IDs, populations and abbreviations
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#sample	pop
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

# Ensure the columns used for matching are correctly named
names(popmap) <- c("individual", "population") # Standardize names if they differ from the file layout

# Assign population factor to genind "pop" based on metadata (outlier name)
pop(outlier.genind)<- as.factor(popmap$population)
pop(outlier.genind)

## levels(outlier.genind$pop) <- c("K", "PP", "CC", "M", "S", "LJP", "TS", "SR") #Based on order in popmap_ref_92

poppr_geneclone <- as.genclone(outlier.genind)

# Finding private alleles by population
private_alleles_result <- poppr::private_alleles(poppr_geneclone, level = "population")

# Convert the result to a data frame
private_alleles_df <- as.data.frame(private_alleles_result)
View(private_alleles_df)

#Don't run (print)
# Function to display private alleles by population
display_private_alleles <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    cat("\nPrivate alleles for population:", pop, "\n")
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print the private alleles and their count
    if (length(private_alleles_pop) > 0) {
      print(private_alleles_pop)
      cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
    } else {
      cat("No private alleles for population", pop, "\n")
    }
  }
}

#Don't run (print)
# Assuming private_alleles_result is your input data (e.g., allele frequencies or counts)
# Example of how to call the function
display_private_alleles(private_alleles_result)

#Run (counts)
# Function to display only the count of private alleles by population
display_private_allele_counts <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print only the count of private alleles
    cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
  }
}

#Run (counts)
# Example of how to call the function
display_private_allele_counts(private_alleles_result)

#Number of private alleles for population Cape Cross : 1 
#Number of private alleles for population False Bay : 1 
#Number of private alleles for population Kleinzee : 1 
#Number of private alleles for population Lambert's Bay : 1 
#Number of private alleles for population Pelican Point : 0 


