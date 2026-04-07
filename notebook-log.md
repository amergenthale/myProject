---
title: "Phylogenetic Analysis of Conserved Genes in Kerstersia gyiorum"
output: html_notebook
---

# Objective

This notebook analyzes phylogenetic relationships among five *Kerstersia gyiorum* genome assemblies using three conserved genes (*rpsJ*, *rplC*, and *rplD*). The workflow includes sequence preparation, multiple sequence alignment, phylogenetic inference, and Bayesian analysis using MrBayes.

# Dataset Description

This project uses a publicly available dataset from a study on bacteria associated with sloths.

**Source:** https://peerj.com/articles/17206/

The organism analyzed is *Kerstersia gyiorum*, a bacterial species isolated from sloths. This dataset contains multiple genome assemblies of the same species.

A subset of five genome assemblies was selected:

- JALJYH000000000  
- JALJYL000000000  
- JALJYN000000000  
- JALJYO000000000  
- JAOQMY000000000  

Three conserved genes were selected for phylogenetic analysis:

- rpsJ  
- rplC  
- rplD  

These genes are commonly used in bacterial phylogenetics and are present across all selected genomes.

# Data Organization

```text
data/
├── rplC.fasta
├── rpsJ.fasta
├── rplD.fasta
```

Data Preparation

The FASTA files were renamed to ensure proper formatting:
```text
mv rplC.fasta.txt rplC.fasta
mv rpsJ.fasta.txt rpsJ.fasta
mv rplD.fasta.txt rplD.fasta
```

Multiple Sequence Alignment

Three alignment methods were tested: ClustalW, MUSCLE, and MAFFT.

Clustal W

clustalw2 -ALIGN -INFILE=rplC.fasta -OUTFILE=rplC-aligned.fasta -OUTPUT=FASTA

Description:
ClustalW is a progressive alignment method that builds a guide tree and aligns sequences stepwise.
Assumptions:

Early alignment decisions are correct
Sequence similarity reflects evolutionary relationships

Limitations:

Errors early in alignment cannot be corrected
Sensitive to guide tree accuracy

MUSCLE
muscle -align rplC.fasta -output rplC-muscle.fasta

escription:
MUSCLE uses iterative refinement and progressive alignment.

Assumptions:

Similar sequences align consistently across iterations

Limitations:

May produce slightly different alignments depending on parameters
Computational cost increases with refinement

MAFFT

mafft --auto rplC.fasta > rplC-mafft.fasta

Description:
MAFFT uses fast Fourier transform and iterative refinement.

Assumptions:

Sequence similarity reflects evolutionary homology
Gap penalties model insertions/deletions

Limitations:

Alignment varies by algorithm mode
Sensitive to highly divergent regions
Alignment Comparison

Three alignment methods (ClustalW, MUSCLE, and MAFFT) were applied to the rplC gene sequences. All methods produced highly similar alignments across conserved regions, indicating strong sequence similarity among taxa. Minor differences occurred in regions with insertions and deletions, where MAFFT and MUSCLE showed different gap placements compared to ClustalW. These differences reflect sensitivity to gap penalties. Regions with inconsistent alignment were treated cautiously. MAFFT was selected for downstream analysis due to its balance of speed and accuracy.

Phylogenetic Inference

All analyses used the MAFFT alignment.

Neighbor-Joining
library(ape)

dna <- read.dna("data/rplC-mafft.fasta", format = "fasta")
dist_matrix <- dist.dna(dna, model = "JC69")

nj_tree <- nj(dist_matrix)

plot(nj_tree, main = "Neighbor-Joining Tree")
Description:
Constructs a tree from pairwise genetic distances.

Assumptions:

Distances reflect evolutionary divergence

Limitations:

Loss of information from sequence simplification
Sensitive to distance errors
Maximum Parsimony
library(phangorn)

dna_phydat <- phyDat(dna, type = "DNA")

start_tree <- nj(dist_matrix)

mp_tree <- optim.parsimony(start_tree, dna_phydat)

plot(mp_tree, main = "Maximum Parsimony Tree")
Description:
Identifies the tree minimizing total evolutionary changes.

Assumptions:

Evolution follows the simplest path

Limitations:

Sensitive to long-branch attraction
No explicit evolutionary model
Maximum Likelihood (IQ-TREE)
cd ~/Desktop/iqtree-3.1.0-macOS/bin
./iqtree3 -s ~/Desktop/myProject/data/rplC-mafft.fasta -bb 1000
Description:
Estimates the most likely tree under a model of evolution.

Assumptions:

Sites evolve independently
Model accurately reflects evolution
Alignment is correct

Limitations:

Computationally intensive
Depends on model selection
Can get trapped in local optima

Tree Visualization and Rooting
library(ape)

tre <- read.tree("data/rplC-mafft.fasta.treefile")

plot(tre)
nodelabels()

rtre <- root(tre, outgroup = "JALJYH000000000", resolve.root = TRUE)

plot(rtre)
nodelabels(rtre$node.label)

Tree Comparison and Interpretation

Neighbor-Joining, Maximum Parsimony, and Maximum Likelihood trees showed broadly similar topologies, indicating consistent evolutionary relationships. Minor differences were observed in branching order, especially in regions with lower support. Maximum Likelihood was considered the most reliable due to its statistical framework, while parsimony may be affected by long-branch attraction.



