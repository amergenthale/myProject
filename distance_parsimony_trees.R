############################################################
# Distance-based and Parsimony-based Phylogenetic Trees
# Author: Audrey Mergenthaler
# Date: 
#
# Description:
# This script estimates:
# 1) A distance-based tree using Neighbor-Joining (NJ)
# 2) A parsimony-based tree using Maximum Parsony (MP)
#
# Packages used:
# - ape
# - phangorn
############################################################

# -----------------------------
# Load Required Packages
# -----------------------------
library(ape)
library(phangorn)

# -----------------------------
# Load Data
# -----------------------------
# Replace with your alignment file
dna <- read.dna("your_alignment.fasta", format = "fasta")

# Convert to phyDat format for phangorn
dna_phydat <- phyDat(dna, type = "DNA")

############################################################
# PART 1: Distance-Based Tree (Neighbor-Joining)
############################################################

# Algorithm Description:
# Neighbor-Joining (NJ) is a distance-based method that
# builds a tree from a matrix of pairwise genetic distances.
# It is fast and computationally efficient.

# Assumptions:
# - Distances accurately represent evolutionary divergence
# - Evolution is approximately additive
# - No explicit model of site-specific evolution

# Limitations:
# - Reduces sequence data to pairwise distances (information loss)
# - Can be inaccurate if model assumptions are violated
# - Sensitive to long-branch attraction

# Calculate distance matrix (Jukes-Cantor model)
dist_matrix <- dist.dna(dna, model = "JC69")

# Build NJ tree
nj_tree <- nj(dist_matrix)

# Plot tree
plot(nj_tree, main="Neighbor-Joining Tree")


############################################################
# PART 2: Parsimony-Based Tree (Maximum Parsimony)
############################################################

# Algorithm Description:
# Maximum Parsimony (MP) finds the tree that minimizes
# the total number of evolutionary changes.

# Assumptions:
# - Evolution follows the simplest path (fewest changes)
# - All substitutions are weighted equally (unless specified)
# - No explicit probabilistic model

# Limitations:
# - Sensitive to long-branch attraction
# - Computationally intensive for large datasets
# - Does not model rate variation among sites

# Create starting tree using NJ
start_tree <- nj(dist_matrix)

# Optimize tree under parsimony
mp_tree <- optim.parsimony(start_tree, dna_phydat)

# Plot MP tree
plot(mp_tree, main="Maximum Parsony Tree")

############################################################
# End of Script
############################################################