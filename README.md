# MQ-PPQ tolerance

This repository includes the scripts developed for the study: "Tolerance of <i>Plasmodium falciparum</i> mefloquine-resistant clinical isolates to paired mefloquine-piperaquine: implications for triple artemisinin-based combination therapy strategies.

# Prerequisites
- Linux Ubuntu 20.04 LTS (or higher)
- Python 3.10.12 (or higher)
- R 4.1 (or higher)

# Producing pileup files from sorted.bam files
Once the reads were mapped to the reference genome, we used the following commands to produce pileup files (one specifically for SNPs, one specifically for indels):

```
#for SNPs
XXX
#for indels
XXX
```

# Counting the A, T, G, C, insertions and deletions along the genome
We then used python scripts to count the different nucleotides, insertions and deletions along the genome:

```
#for SNPs
XXX
#for indels
XXX
```

# Comparison of the two genomes
We finally developped two R scripts to search for mutations (SNPs, indels) in the paired samples. Some filtration steps are included to remove false positive variants, and not consider variants that could be present at minor proportions in one sample and at major proportions in the other one. The script must be executed in Rstudio.

# Citation
**XXX**

Camille Roesh, Romain Coppée, etc.
