# Reproducible Notebook

## Dataset Description
I am using a publicly available dataset: https://peerj.com/articles/17206/

This dataset contains analysis of twelve genomes of the bacterium Kerstersia gyiorum from brown-throated sloths (Bradypus variegatus), the first from a non-human host

The dataset was selected because it is relevant to evolutionary biology and allows downstream phylogenetic analyses without requiring raw sequencing data.

## Data Source
This dataset conatins 12 species, I am choosing 5 species and 3 genes from each. I will organize the data into Gene rpsJ, rplC, rplD. 

- Repository: myProject
- Organism(s): Sloths 

## Steps Taken
- Downloaded 5 species data sets based on the 3 genes I chose to anaylyze 
- The 5 species are: JALJYH000000000, JALJYL000000000, JALJYN000000000, JALJYO000000000, JAOQMY000000000
- The 3 genes I chose to analyze are: rpsJ, rplC, rplD 
- Organized each gene into its own file  

## Clustal W
COMMAND USED: clustalw2 -ALIGN -INFILE=rplC.fasta -OUTFILE=rplC-aligned.fasta -OUTPUT=FASTA

CLUSTAL 2.1 Multiple Sequence Alignments


Sequence format is Pearson
Sequence 1: NZ_JALJYH010000001.1_3687-4376       690 bp
Sequence 2: NZ_JALJYL010000001.1_656747-657436   690 bp
Sequence 3: NZ_JALJYN010000001.1_656733-657422   690 bp
Sequence 4: NZ_JALJYO010000001.1_656740-657429   690 bp
Sequence 5: NZ_JAOQMY010000001.1_656733-657422   690 bp
Start of Pairwise alignments
Aligning...

Sequences (1:2) Aligned. Score:  11
Sequences (1:3) Aligned. Score:  11
Sequences (1:4) Aligned. Score:  11
Sequences (1:5) Aligned. Score:  11
Sequences (2:3) Aligned. Score:  100
Sequences (2:4) Aligned. Score:  99
Sequences (2:5) Aligned. Score:  100
Sequences (3:4) Aligned. Score:  99
Sequences (3:5) Aligned. Score:  100
Sequences (4:5) Aligned. Score:  99
Guide tree file created:   [rplC.fasta.dnd]

There are 4 groups
Start of Multiple Alignment

Aligning...
Group 1: Sequences:   2      Score:13110
Group 2: Sequences:   3      Score:13110
Group 3: Sequences:   4      Score:13091
Group 4:                     Delayed
Alignment Score 40090
firstres = 1 lastres = 718
FASTA file created!

Fasta-Alignment file created    [rplC-aligned.fasta]



