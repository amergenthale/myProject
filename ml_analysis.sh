#!/bin/bash

# ==========================================
# Maximum Likelihood Phylogenetic Analysis
# Using IQ-TREE
# ==========================================

# DESCRIPTION:
# IQ-TREE is a Maximum Likelihood (ML) method that estimates
# phylogenetic trees by finding the tree topology that maximizes
# the likelihood of observing the given sequence data.

# ASSUMPTIONS:
# - Sites evolve independently
# - A substitution model (e.g., GTR) describes sequence evolution
# - The sequence alignment is correct

# LIMITATIONS:
# - Computationally intensive
# - Can get stuck in local optima
# - Results depend on model selection

# COMMAND (not executed for this assignment):
iqtree -s alignment.fasta -m TEST -bb 1000 -nt AUTO