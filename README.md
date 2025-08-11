# MQ-PPQ tolerance

This repository includes the scripts developed for the study: "Tolerance of <i>Plasmodium falciparum</i> mefloquine-resistant clinical isolates to paired mefloquine-piperaquine: implications for triple artemisinin-based combination therapy strategies.

# Prerequisites
- Linux Ubuntu 20.04 LTS (or higher)
- Python 3.10.12 (or higher)
- R 4.1 (or higher)

# Producing pileup files from sorted.bam files
Once the reads were mapped to the reference genome, we used the following commands to produce pileup files (one specifically for SNPs, one specifically for indels):

```
##for SNPs, example for sample F75
samtools mpileup -a -A -x -B --no-output-ins --no-output-del -f Pfalciparum.genome.fasta F75.sorted.bam -o F75-all.pileup
##for indels, remove --no-output-ins and --no-output-del parameters
```

# Counting the A, T, G, C, insertions and deletions along the genome
We then used python scripts to count the different nucleotides, insertions and deletions along the genome:

```
##for SNPs (change the file name in the script)
python3 pileup_to_count_bases.py
##for indels
python3 pileup_to_count_indels.py
```

# Comparison of the two genomes
We finally developped two R scripts to search for mutations (SNPs or indels) in the paired samples. Some filtration steps are included to remove false positive variants, and not consider variants that could be present at minor proportions in one sample and at major proportions in the other one. The script must be executed in Rstudio.

# Citation
**to do**

Camille Roesh, Romain Coppée, TO COMPLETE etc.
