#!/bin/bash
#SBATCH --job-name="Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s PB2_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s PB1_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s PA_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s NS_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s NP_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s NA_N2_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s NA_N1_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s MP_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s HA_H3_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s HA_H1_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
