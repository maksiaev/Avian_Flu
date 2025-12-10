#!/bin/bash
#SBATCH --job-name="all random iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s 1-random_PB2_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 2-random_PB1_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 3-random_PA_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 4-random_HA_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 5-random_NP_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 6-random_NA_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 7-random_MP_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO
iqtree2 -s 8-random_NS_2021-11-01--2025-12-05_aln.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
