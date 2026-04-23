---
Title: "Phylogenetic Analysis of Conserved Genes in Kerstersia gyiorum"
---

# Objective

his notebook aims to investigate the evolutionary relationships among five *Kerstersia gyiorum* genome assemblies using two conserved genes (*rpsJ* and *rplC*). By analyzing multiple genes, this study provides a more robust assessment of phylogenetic patterns within the species. 

# Dataset Description

This project uses a publicly available dataset from a study on bacteria associated with sloths.

**Source:** https://peerj.com/articles/17206/

The organism analyzed is *Kerstersia gyiorum*, a Gram-negative bacterium found in both humans and animals and associated with opportunistic infections. The genomes used in this analysis were isolated from brown-throated sloths (*Bradypus variegatus*), an arboreal mammal native to Central and South America. Sloths have a unique lifestyle and slow metabolism, which contribute to a distinct microbiome. Studying bacteria from sloths provides insight into host-associated microbial diversity and potential evolutionary adaptation to different environments. 

Although the original study analyzed multiple genomes, a subset of five genome assemblies was selected for this analysis:

A subset of five genome assemblies was selected:

- JALJYH000000000  
- JALJYL000000000  
- JALJYN000000000  
- JALJYO000000000  
- JAOQMY000000000  

Three conserved genes were selected for phylogenetic analysis:

- rpsJ (ribosomal protein S10)  
- rplC (ribosomal protein L3)   

These genes encode components of the bacterial ribosome and are essential for protein synthesis. Because they are highly conserved across bacteria but still accumulate mutations over time, they are well-suited for reconstructing evolutionary relationships. Their presence in all genomes ensures comparability across taxa, while their variation provides phylogenetic signal.

# Data Organization

```text
data/
├── rplC.fasta
├── rpsJ.fasta
```

# Data Preparation

The FASTA files were renamed to ensure proper formatting:
```text
mv rplC.fasta.txt rplC.fasta
mv rpsJ.fasta.txt rpsJ.fasta
```

# Multiple Sequence Alignment

Three alignment methods were tested: ClustalW, MUSCLE, and MAFFT.

# Clustal W

```text
clustalw2 -ALIGN -INFILE=rplC.fasta -OUTFILE=rplC-aligned.fasta -OUTPUT=FASTA
```

**Description:**
ClustalW is a progressive alignment method that builds a guide tree and aligns sequences stepwise.

**Assumptions:**
Early alignment decisions are correct
Sequence similarity reflects evolutionary relationships

**Limitations:**
Errors early in alignment cannot be corrected
Sensitive to guide tree accuracy

# MUSCLE
```text
muscle -align rplC.fasta -output rplC-muscle.fasta
```

**Description:**
MUSCLE uses iterative refinement and progressive alignment.

**Assumptions:**
Similar sequences align consistently across iterations

**Limitations:**
May produce slightly different alignments depending on parameters
Computational cost increases with refinement

# MAFFT
```text
mafft --auto rplC.fasta > rplC-mafft.fasta
```

**Description:**
MAFFT uses fast Fourier transform and iterative refinement.

**Assumptions:**
Sequence similarity reflects evolutionary homology
Gap penalties model insertions/deletions

**Limitations:**
Alignment varies by algorithm mode
Sensitive to highly divergent regions

# Alignment Comparison
Three alignment methods (ClustalW, MUSCLE, and MAFFT) were applied to the rplC gene sequences. All methods produced highly similar alignments across conserved regions, indicating strong sequence similarity among taxa. Minor differences occurred in regions with insertions and deletions, where MAFFT and MUSCLE showed different gap placements compared to ClustalW. These differences reflect sensitivity to gap penalties. Regions with inconsistent alignment were treated cautiously. MAFFT was selected for downstream analysis due to its balance of speed and accuracy.

# Phylogenetic Inference

All analyses used the MAFFT alignment.

# Neighbor-Joining
```text
library(ape)
dna <- read.dna("data/rplC-mafft.fasta", format = "fasta")
dist_matrix <- dist.dna(dna, model = "JC69")
nj_tree <- nj(dist_matrix)
plot(nj_tree, main = "Neighbor-Joining Tree")
```
**Description:**
Constructs a tree from pairwise genetic distances.

**Assumptions:**
Distances reflect evolutionary divergence

**Limitations**:
Loss of information from sequence simplification
Sensitive to distance errors

# Maximum Parsimony
```text
library(phangorn)
dna_phydat <- phyDat(dna, type = "DNA")
start_tree <- nj(dist_matrix)
mp_tree <- optim.parsimony(start_tree, dna_phydat)
plot(mp_tree, main = "Maximum Parsimony Tree")
```

**Description:**
Identifies the tree minimizing total evolutionary changes.

**Assumptions:**
Evolution follows the simplest path

**Limitations:**
Sensitive to long-branch attraction
No explicit evolutionary model

# Maximum Likelihood (IQ-TREE)
```text
cd ~/Desktop/iqtree-3.1.0-macOS/bin
./iqtree3 -s ~/Desktop/myProject/data/rplC-mafft.fasta -bb 1000
```

**Description:**
Estimates the most likely tree under a model of evolution.

**Assumptions:**
Sites evolve independently
Model accurately reflects evolution
Alignment is correct

**Limitations:**
Computationally intensive
Depends on model selection
Can get trapped in local optima

# Tree Visualization and Rooting
```text
library(ape)
tre <- read.tree("data/rplC-mafft.fasta.treefile")
plot(tre)
nodelabels()
rtre <- root(tre, outgroup = "JALJYH000000000", resolve.root = TRUE)
plot(rtre)
nodelabels(rtre$node.label)
```

# Tree Comparison and Interpretation
Neighbor-Joining, Maximum Parsimony, and Maximum Likelihood trees showed broadly similar topologies, indicating consistent evolutionary relationships. Minor differences were observed in branching order, especially in regions with lower support. Maximum Likelihood was considered the most reliable due to its statistical framework, while parsimony may be affected by long-branch attraction.

# Bayesian Phylogenetic Inference with MrBayes

**Description:**
Bayesian phylogenetic inference was performed using MrBayes. MrBayes uses a Markov Chain Monte Carlo (MCMC) algorithm to estimate the posterior probability distribution of trees under a specified model of sequence evolution. Instead of producing a single tree, it samples many trees according to their probability given the data and model.

**Assumptions:**
The chosen substitution model (GTR + Γ) adequately represents sequence evolution
Sites evolve independently
The alignment is correct
The Markov chain reaches convergence and samples the true posterior distribution

**Limitations:**
Results depend heavily on the chosen model
Poor alignments can lead to incorrect trees
MCMC may not converge if run too briefly
Computationally intensive for large datasets

MrBayes Commands
```text
mb
execute rplC-mafft-clean.nex
lset nst=6 rates=gamma
mcmc ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 nchains=4
sump burnin=2500
sumt burnin=2500
```

**Notes:**
Sequence alignment was performed using MAFFT
FASTA files were converted to NEXUS format using the R package ape
Taxon names were cleaned to remove spaces and special characters for compatibility with MrBayes
Errors encountered during file formatting were corrected (e.g., missing ;, missing end;, invalid taxon names)

# Coalescent Species Tree Analysis
**Goal**
Because this dataset contains three genes (`rpsJ` and `rplC`), a coalescent-based species tree method can be applied to the project data. Unlike single-gene tree inference, coalescent methods use multiple gene trees to estimate a species tree while accounting for discordance among loci.

**Chosen Method**
I chose **ASTRAL** as the coalescent method.

**Description of the Algorithm
ASTRAL is a summary coalescent method that estimates a species tree from a set of input gene trees. Rather than analyzing raw sequence alignments directly, it uses the topologies of gene trees and identifies the species tree that is most consistent with the quartet relationships observed across loci. This makes it useful when different genes support somewhat different evolutionary histories.

**Assumptions**
- Multiple independent loci are available.
- Gene tree discordance is mainly caused by incomplete lineage sorting.
- Input gene trees are reasonably accurate.
- Genes are orthologous rather than paralogous.
- Recombination within each locus is low, and loci are treated as independent.

**Limitations**
- ASTRAL depends on the quality of the input gene trees.
- With only a small number of genes, species tree estimation may have limited resolution.
- It is designed mainly for discordance caused by incomplete lineage sorting and does not explicitly model processes such as horizontal gene transfer or hybridization.
- Errors in alignment or gene tree estimation can affect the final species tree.

**Align each gene**
```bash
mafft --auto data/rpsJ.fasta > data/rpsJ-mafft.fasta
mafft --auto data/rplC.fasta > data/rplC-mafft.fasta
```

**Infer ML gene trees**
```bash
iqtree2 -s data/rpsJ-mafft.fasta -m MFP -bb 1000 -nt AUTO
iqtree2 -s data/rplC-mafft.fasta -m MFP -bb 1000 -nt AUTO
```

**Combine gene trees**
```bash
cat data/rpsJ-mafft.fasta.treefile \
    data/rplC-mafft.fasta.treefile \
```

**Run ASTRAL**
```bash
java -jar astral.jar -i data/all_gene_trees.tre -o data/astral_species_tree.tre
```

**Notes** 
This coalescent workflow uses the three genes in the project dataset rather than a toy example. The assignment only requires the commands and method description, so obtaining the final species tree output is not necessary for submission.









