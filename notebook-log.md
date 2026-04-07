Reproducible Notebook 

DATA DESCRIPTION
This project uses a publicly available dataset from a study on the bacterium associated with sloths:
Source: https://peerj.com/articles/17206/
Organism: Kerstersia gyiorum (bacterial genomes isolated from sloths)
This dataset contains multiple genome assemblies of the same bacterial species. A subset of five genomes was selected for analysis:
JALJYH000000000
JALJYL000000000
JALJYN000000000
JALJYO000000000
JAOQMY000000000
Three conserved genes were selected for phylogenetic analysis:
rpsJ
rplC
rplD
These genes were chosen because they are commonly used in bacterial phylogenetics and are present across all selected genomes.

DATA ORGANIZACTION
Each gene was extracted and stored as a separate FASTA file:
data/
├── rplC.fasta
├── rpsJ.fasta
├── rplD.fasta

DATA PREPARATION
One issue encountered during preprocessing was incorrect file extensions (.fasta.txt). These were corrected:
mv rplC.fasta.txt rplC.fasta
mv rpsJ.fasta.txt rpsJ.fasta
mv rplD.fasta.txt rplD.fasta
This step was necessary because alignment software requires proper FASTA formatting and file naming.
Multiple Sequence Alignment
Three alignment methods were used to compare performance and ensure robustness:

CLUSTAL W
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
Description:
MUSCLE uses iterative refinement and progressive alignment to improve accuracy.
Assumptions:
Similar sequences should align consistently across iterations
Limitations:
May produce slightly different alignments depending on parameters
Computational cost increases with refinement

MAFT
mafft --auto rplC.fasta > rplC-mafft.fasta
Description:
MAFFT uses fast Fourier transform techniques and iterative refinement to efficiently align sequences.
Assumptions:
Sequence similarity reflects evolutionary homology
Gap penalties appropriately model insertions/deletions
Limitations:
Alignment can vary depending on chosen algorithm mode
Sensitive to highly divergent regions

ALIGNMENT COMPARISON
Three alignment methods (ClustalW, MUSCLE, and MAFFT) were applied to the rplC gene sequences. Overall, all three methods produced highly similar alignments across conserved regions, indicating strong sequence similarity among the selected taxa. However, minor differences were observed in regions containing insertions and deletions, where MAFFT and MUSCLE tended to introduce slightly different gap placements compared to ClustalW. These differences highlight the sensitivity of alignment algorithms to gap penalties and refinement strategies. Regions with inconsistent gap placement across methods were considered unreliable and were interpreted cautiously in downstream phylogenetic analyses. Because MAFFT balances speed and accuracy and is widely recommended for nucleotide data, its alignment was selected for subsequent tree inference.

PHLYOGENETIC INFERENCE
Distance-Based Method (Neighbor-Joining)
library(ape)
dna <- read.dna("rplC-mafft.fasta", format="fasta")
dist_matrix <- dist.dna(dna, model="JC69")
nj_tree <- nj(dist_matrix)
plot(nj_tree, main="Neighbor-Joining Tree")
Description:
Neighbor-Joining constructs a tree from pairwise distances.
Assumptions:
Distances accurately reflect evolutionary divergence
Limitations:
Loss of information due to distance conversion
Sensitive to inaccurate distance estimation

PARISMONY METHOD (MAXIMUM PARSIMONY)
library(phangorn)
dna_phydat <- phyDat(dna, type="DNA")
start_tree <- nj(dist_matrix)
mp_tree <- optim.parsimony(start_tree, dna_phydat)
plot(mp_tree, main="Maximum Parsimony Tree")
Description:
Maximum Parsimony identifies the tree requiring the fewest evolutionary changes.
Assumptions:
Evolution follows the simplest path
Limitations:
Sensitive to long-branch attraction
Does not use an explicit evolutionary model

MAXIMUM LIKELIHOOD (IG-TREE)
cd ~/Desktop/iqtree-3.1.0-macOS/bin
./iqtree3 -s ~/Desktop/myProject/data/rplC-mafft.fasta -bb 1000
Description:
Maximum Likelihood estimates the tree that maximizes the probability of the observed data under a model of evolution.
Assumptions:
Sites evolve independently
The substitution model accurately reflects sequence evolution
Alignment is correct
Limitations:
Computationally intensive
Can be sensitive to model choice
Plot ML Tree
library(ape)
tre <- read.tree("rplC-mafft.fasta.treefile")
plot(tre)
nodelabels()
Rooting the Tree
rtre <- root(tre, outgroup="JALJYH000000000", resolve.root=TRUE)
plot(rtre)
nodelabels (rtre$node.label)

TREE COMPARISON
The Neighbor-Joining, Maximum Parsimony, and Maximum Likelihood trees showed broadly similar topologies, indicating consistent evolutionary relationships among the sampled taxa. However, minor differences in branching order were observed, particularly in regions with lower support values. The Maximum Likelihood tree was considered the most reliable due to its use of an explicit evolutionary model and statistical framework, while the parsimony method may be affected by long-branch attraction. These differences highlight the importance of comparing multiple phylogenetic methods.




