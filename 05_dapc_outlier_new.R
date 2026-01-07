############################################################
# 0. HOUSEKEEPING
############################################################

# House keeping ####
#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts")

# Generate citations ####
#For R
citation()

#For R Studio
if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")
library(rstudioapi)

versionInfo()$citation

#For packages
#citation("pkgname")
# Load packages ####
if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
if (!require("vcfR", quietly = TRUE)) install.packages("vcfR")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tibble", quietly = TRUE)) install.packages("tibble")
if (!require("scales", quietly = TRUE)) install.packages("scales")
if (!require("viridis", quietly = TRUE)) install.packages("viridis")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(adegenet)
library(vcfR)
library(dplyr)
library(tibble)
library(scales)
library(viridis)
library(ggplot2)

# 1. READ DATA
outlier.vcf <- read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_outlier.vcf")
outlier.genlight <- vcfR2genlight(outlier.vcf)

# 2. FIND CLUSTERS (INTERACTIVE)
num_clust.outlier <- find.clusters(
  outlier.genlight, 
  scale = FALSE
)
# number of PCs to retain (can be all, see manual) 55 (I used the no of PCs explaining aprox ~85% variance) 
# number of clusters (lowest BIC)

# OR if you have alread done the find.cluster step: 
#num_clust.outlier <- find.clusters(
#  outlier.genlight, 
# scale = FALSE,
#  n.pca = 55,       # Initial PCs for K-means
#  n.clust = 2       # Optimal K from BIC plot
#)

# 3. ASSIGN GROUPS TO GENLIGHT OBJECT (CRITICAL STEP)
# Assign the K-means groups (K=2) to the population slot of the genlight object
# NOTE: This overwrites the field population set in Step 1c, which is fine for DAPC grouping.
pop(outlier.genlight) <- num_clust.outlier$grp

# 4. INITIAL DAPC FOR OPTIMIZATION
# Rerun DAPC on the now-ordered genlight object
outlier.dapc_initial <- dapc(
  outlier.genlight, 
  n.pca = 55, 
  n.da = 1 
)

# 5. OPTIMIZE n.pca 
opt.a.score <- optim.a.score(outlier.dapc_initial, plot = FALSE)
best_n_pca <- opt.a.score$best 
best_n_pca
# Result is 1 
# Agrees with Jombart and Collins (2023), which says this should NOT be > K-1 

# 6. FINAL DAPC RUN
# Rerun FINAL DAPC on the now-ordered genlight object
outlier.dapc_final <- dapc(
  outlier.genlight, 
  n.pca = best_n_pca, 
  n.da = 1 
)



# Load your popmap for plotting only
popmap<-read.delim("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/popmap91", header = TRUE, sep = "\t")
#file layout:
#individual	population
#CC01 Cape Cross
#CC03 Cape Cross
#CC04 Cape Cross
#...
names(popmap) <- c("individual", "population")
popmap$individual <- gsub("^KZ", "KL", popmap$individual)

# Keep only individuals present in the genlight
popmap_filtered <- popmap %>% filter(individual %in% indNames(outlier.genlight))

# Define plotting order
pop_order_custom <- c("Cape Cross", "Pelican Point", "Kleinzee", "Lambert's Bay", "False Bay")
popmap_filtered <- popmap_filtered %>%
  mutate(population = factor(population, levels = pop_order_custom)) %>%
  arrange(population)

# Add ordered colony info to 'other' slot for plotting
other(outlier.genlight)$pop <- popmap_filtered$population[match(indNames(outlier.genlight), popmap_filtered$individual)]

# 8. DEFINE COLOR PALETTE 
# IMPORTANT: The number of colors must match the number of UNIQUE populations (5 in this case)
palette_5pop <- c(
  "Cape Cross" = "#414487FF",
  "Pelican Point" = "#440154FF", # Moved Pelican Point up
  "Kleinzee" = "#2A788EFF",
  "Lambert's Bay" = "#22A884FF",
  "False Bay" = "#7AD151FF" # Moved False Bay down
) 

# 1️⃣ Define your plotting colony order (same as old code)
colony_order <- c("Cape Cross", "Pelican Point", "Kleinzee", "Lambert's Bay", "False Bay")

# 2️⃣ Ensure 'other' slot is populated with colonies
other(outlier.genlight)$pop <- factor(
  popmap_filtered$population[match(indNames(outlier.genlight), popmap_filtered$individual)],
  levels = colony_order  # enforce correct plotting order
)

# 3️⃣ Prepare DAPC coordinates for plotting
dapc_coords_df <- as_tibble(outlier.dapc_final$ind.coord, rownames = "individual") %>%
  mutate(
    group = factor(outlier.dapc_final$grp),      # de novo DAPC groups
    population = other(outlier.genlight)$pop     # colony for plotting
  )

# 4️⃣ 1D Density Plot (fill = DAPC groups, not colony)
ggplot(dapc_coords_df, aes(x = LD1, fill = group)) +
  geom_density(alpha = 0.6) +
  labs(
    title = "DAPC 1D Density Plot (K=2 Outliers)",
    subtitle = "Individuals colored by DAPC-inferred group",
    x = "Linear Discriminant 1 (LD1)",
    y = "Density"
  ) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red")) +
  theme_minimal(base_size = 14)




# Define desired colony order
pop_order_custom <- c("Cape Cross", "Pelican Point", "Kleinzee", "Lambert's Bay", "False Bay")

# Filter popmap to include only individuals present in genlight
popmap_filtered <- popmap %>%
  filter(individual %in% indNames(outlier.genlight)) %>%
  mutate(population = factor(population, levels = pop_order_custom)) %>%
  arrange(population)

# Reorder genlight numerically
match_indices <- match(popmap_filtered$individual, indNames(outlier.genlight))
outlier.genlight <- outlier.genlight[match_indices]

# Assign population factor to 'pop' slot (optional: used for plotting)
pop(outlier.genlight) <- popmap_filtered$population

# Store original colony info in 'other' slot for compoplot labeling
other(outlier.genlight)$pop <- popmap_filtered$population

# --- RUN FIND.CLUSTERS AND DAPC --- #
num_clust.outlier <- find.clusters(outlier.genlight, scale = FALSE, n.pca = 55, n.clust = 2)
pop(outlier.genlight) <- num_clust.outlier$grp  # assign DAPC groups

outlier.dapc_final <- dapc(outlier.genlight, n.pca = 1, n.da = 1)  # using optimized n.pca

# --- COMPOPLOT WITH REORDERED BARS --- #
compoplot(
  outlier.dapc_final,
  col = c("blue", "red"),                   # DAPC-inferred groups
  lab = indNames(outlier.genlight),
  show.lab = TRUE,
  cex.names = 0.3,
  group.names = other(outlier.genlight)$pop, # ORIGINAL colony labels in desired order
  grp.leg = FALSE,
  spacing = 0.15,
  xlab = "Individuals Grouped by Field Population (Reordered)",
  legend = FALSE,
  cex.lab = 1.2
)


############################################################
# DETERMINE COLORS TO MATCH COMPOPLOT
############################################################

# Extract cluster assignments in the plotting order
cluster_ordered <- dapc_df$cluster

# Create color vector: use same colors as compoplot
cluster_colors <- c("1" = "blue", "2" = "red")  # swap if needed to match compoplot

############################################################
# DENSITY PLOT MATCHING COMPOPLOT COLORS
############################################################

p_density <- ggplot(dapc_df, aes(x = LD1, fill = cluster)) +
  geom_density(alpha = 0.6) +
  labs(
    title = "DAPC Density Plot — Clusters Matching Compoplot",
    x = "LD1",
    y = "Density"
  ) +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal(base_size = 14)

print(p_density)

############################################################
# COMPOPLOT WITH INDIVIDUALS ORDERED BY FIELD POPULATION
############################################################

compoplot(
  full.dapc,
  col = c("blue", "red"),  # same colors as density
  lab = dapc_df$individual,
  show.lab = FALSE,
  legend = TRUE,
  cex.names = 0.5,
  grp.leg = FALSE,
  spacing = 0.15  # space between populations
)